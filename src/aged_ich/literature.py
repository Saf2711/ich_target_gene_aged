from __future__ import annotations

import csv
import gzip
import json
import re
import time
import xml.etree.ElementTree as ET
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterator, Sequence

import pandas as pd
import requests
from tqdm import tqdm

from .config import (
    DEFAULT_ICH_QUERY,
    EUTILS_BASE,
    NCBI_GENE_INFO_URL,
    PUBTATOR_EXPORT,
)


REVIEW_PUBLICATION_TYPES = {
    "review",
    "systematic review",
    "meta-analysis",
    "editorial",
    "comment",
    "letter",
    "practice guideline",
    "guideline",
}

INCREASE_TERMS = {
    "overexpress",
    "overexpression",
    "upregulat",
    "recombinant",
    "supplement",
    "administration",
    "administered",
    "delivery",
    "delivered",
    "transfection",
    "transduced",
    "gene transfer",
    "viral vector",
    "aav",
    "adenovirus",
    "lentivirus",
    "agonist",
    "activator",
    "activation",
    "extracellular vesicle",
    "exosome",
}

DECREASE_TERMS = {
    "knockout",
    "knock-out",
    "knockdown",
    "knock-down",
    "silencing",
    "silenced",
    "deletion",
    "deficiency",
    "inhibitor",
    "inhibition",
    "inhibited",
    "antagonist",
    "neutralizing antibody",
    "blocked",
    "blockade",
}

BENEFIT_TERMS = {
    "improved",
    "improvement",
    "ameliorated",
    "attenuated",
    "reduced injury",
    "reduced edema",
    "reduced oedema",
    "reduced hematoma",
    "reduced haematoma",
    "reduced neuronal death",
    "neuroprotect",
    "protected",
    "recovery",
    "better neurological",
    "improved neurological",
    "decreased mortality",
    "reduced mortality",
    "enhanced clearance",
    "hematoma clearance",
    "haematoma clearance",
}

HARM_TERMS = {
    "worsened",
    "exacerbated",
    "aggravated",
    "increased injury",
    "increased edema",
    "increased oedema",
    "increased hematoma",
    "increased haematoma",
    "increased neuronal death",
    "neurotoxic",
    "higher mortality",
    "increased mortality",
    "impaired recovery",
}

MECHANISM_TERMS = {
    "mechanism",
    "pathway",
    "mediated",
    "dependent",
    "regulates",
    "regulated",
    "associated with",
    "correlated with",
    "biomarker",
    "expression",
    "protein level",
    "mrna",
    "immunohistochemistry",
    "western blot",
}

MODEL_TERMS = {
    "collagenase",
    "autologous blood",
    "blood injection",
    "hemin",
    "haemin",
    "hemoglobin",
    "haemoglobin",
    "intracerebral hemorrhage model",
    "intracerebral haemorrhage model",
    "mouse",
    "mice",
    "rat",
    "rats",
    "patient",
    "patients",
    "human",
}


@dataclass(frozen=True)
class ApiSettings:
    """NCBI request settings."""

    email: str
    api_key: str | None = None
    tool: str = "aged_ich_target_prioritization"
    timeout_seconds: int = 60

    @property
    def requests_per_second(self) -> float:
        return 10.0 if self.api_key else 3.0

    @property
    def delay_seconds(self) -> float:
        return 1.0 / self.requests_per_second


class NcbiClient:
    """NCBI and PubTator client with retry and rate-limit handling."""

    def __init__(
        self,
        settings: ApiSettings,
        cache_dir: Path,
    ) -> None:
        self.settings = settings
        self.cache_dir = cache_dir
        self.cache_dir.mkdir(parents=True, exist_ok=True)

        self.session = requests.Session()
        self.session.headers.update(
            {
                "User-Agent": (
                    f"{settings.tool}/1.0 ({settings.email})"
                )
            }
        )

        self._last_request_time = 0.0

    def _wait(self) -> None:
        elapsed = time.monotonic() - self._last_request_time
        wait_for = self.settings.delay_seconds - elapsed

        if wait_for > 0:
            time.sleep(wait_for)

    def _request(
        self,
        method: str,
        url: str,
        *,
        params: dict[str, Any] | None = None,
        data: dict[str, Any] | None = None,
        retries: int = 5,
    ) -> requests.Response:
        for attempt in range(retries):
            self._wait()

            try:
                response = self.session.request(
                    method,
                    url,
                    params=params,
                    data=data,
                    timeout=self.settings.timeout_seconds,
                )
                self._last_request_time = time.monotonic()

                if response.status_code in {
                    429,
                    500,
                    502,
                    503,
                    504,
                }:
                    raise requests.HTTPError(
                        f"Temporary HTTP {response.status_code}",
                        response=response,
                    )

                response.raise_for_status()
                return response

            except (
                requests.ConnectionError,
                requests.Timeout,
                requests.HTTPError,
            ):
                if attempt == retries - 1:
                    raise

                time.sleep(min(60.0, 2.0 ** attempt))

        raise RuntimeError("Request retry loop ended unexpectedly.")

    def eutils_params(self) -> dict[str, str]:
        params = {
            "email": self.settings.email,
            "tool": self.settings.tool,
        }

        if self.settings.api_key:
            params["api_key"] = self.settings.api_key

        return params

    def search_pubmed(
        self,
        query: str = DEFAULT_ICH_QUERY,
    ) -> list[str]:
        """Return every PubMed identifier matching the ICH query."""

        count_params = {
            **self.eutils_params(),
            "db": "pubmed",
            "term": query,
            "retmode": "json",
            "retmax": "0",
        }

        count_response = self._request(
            "GET",
            f"{EUTILS_BASE}/esearch.fcgi",
            params=count_params,
        )

        count = int(
            count_response.json()["esearchresult"]["count"]
        )

        pmids: list[str] = []
        page_size = 10000

        for start in tqdm(
            range(0, count, page_size),
            desc="Retrieving ICH PubMed identifiers",
        ):
            params = {
                **self.eutils_params(),
                "db": "pubmed",
                "term": query,
                "retmode": "json",
                "retstart": str(start),
                "retmax": str(page_size),
            }

            response = self._request(
                "GET",
                f"{EUTILS_BASE}/esearch.fcgi",
                params=params,
            )

            pmids.extend(
                response.json()["esearchresult"]["idlist"]
            )

        return list(dict.fromkeys(pmids))

    def fetch_pubmed_xml(
        self,
        pmids: Sequence[str],
    ) -> bytes:
        params = {
            **self.eutils_params(),
            "db": "pubmed",
            "id": ",".join(pmids),
            "retmode": "xml",
            "rettype": "abstract",
        }

        response = self._request(
            "POST",
            f"{EUTILS_BASE}/efetch.fcgi",
            data=params,
        )

        return response.content

    def fetch_pubtator_biocjson(
        self,
        pmids: Sequence[str],
    ) -> Any:
        response = self._request(
            "GET",
            PUBTATOR_EXPORT,
            params={"pmids": ",".join(pmids)},
        )

        return response.json()

    def download(
        self,
        url: str,
        destination: Path,
    ) -> Path:
        if (
            destination.exists()
            and destination.stat().st_size > 0
        ):
            return destination

        destination.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        response = self._request("GET", url)
        destination.write_bytes(response.content)

        return destination


def chunks(
    values: Sequence[str],
    size: int,
) -> Iterator[list[str]]:
    for start in range(0, len(values), size):
        yield list(values[start : start + size])


def text_or_empty(
    element: ET.Element | None,
) -> str:
    if element is None:
        return ""

    return " ".join(
        "".join(element.itertext()).split()
    )


def parse_pubmed_articles(
    xml_content: bytes,
) -> list[dict[str, Any]]:
    """Convert PubMed XML records into a flat article table."""

    root = ET.fromstring(xml_content)
    records: list[dict[str, Any]] = []

    for article_node in root.findall(".//PubmedArticle"):
        citation = article_node.find("./MedlineCitation")

        article = (
            citation.find("./Article")
            if citation is not None
            else None
        )

        if citation is None or article is None:
            continue

        pmid = text_or_empty(
            citation.find("./PMID")
        )

        title = text_or_empty(
            article.find("./ArticleTitle")
        )

        abstract_parts: list[str] = []

        for abstract_text in article.findall(
            "./Abstract/AbstractText"
        ):
            label = abstract_text.attrib.get(
                "Label",
                "",
            )
            text = text_or_empty(abstract_text)

            if label and text:
                abstract_parts.append(
                    f"{label}: {text}"
                )
            elif text:
                abstract_parts.append(text)

        abstract = " ".join(abstract_parts)

        publication_types = [
            text_or_empty(node)
            for node in article.findall(
                "./PublicationTypeList/PublicationType"
            )
            if text_or_empty(node)
        ]

        year = ""

        for path in (
            "./Journal/JournalIssue/PubDate/Year",
            "./ArticleDate/Year",
            "../DateCompleted/Year",
            "../DateRevised/Year",
        ):
            year = text_or_empty(article.find(path))

            if year:
                break

        doi = ""

        for article_id in article_node.findall(
            "./PubmedData/ArticleIdList/ArticleId"
        ):
            if article_id.attrib.get("IdType") == "doi":
                doi = text_or_empty(article_id)
                break

        records.append(
            {
                "pmid": pmid,
                "title": title,
                "abstract": abstract,
                "year": year,
                "publication_types": "|".join(
                    publication_types
                ),
                "doi": doi,
            }
        )

    return records


def iter_bioc_documents(
    payload: Any,
) -> Iterator[dict[str, Any]]:
    """Yield BioC documents from common PubTator response shapes."""

    if isinstance(payload, list):
        for item in payload:
            if (
                isinstance(item, dict)
                and "documents" in item
            ):
                yield from item["documents"]
            elif isinstance(item, dict):
                yield item

        return

    if isinstance(payload, dict):
        if "documents" in payload:
            yield from payload["documents"]

        elif (
            "collection" in payload
            and isinstance(
                payload["collection"],
                dict,
            )
        ):
            yield from payload[
                "collection"
            ].get(
                "documents",
                [],
            )

        else:
            yield payload


def annotation_gene_id(
    annotation: dict[str, Any],
) -> str | None:
    infons = annotation.get("infons", {}) or {}

    gene_id = (
        infons.get("identifier")
        or infons.get("Identifier")
        or infons.get("NCBI Gene")
    )

    if gene_id is None:
        return None

    match = re.search(r"\d+", str(gene_id))
    return match.group(0) if match else None


def parse_pubtator_gene_pairs(
    payload: Any,
) -> list[dict[str, str]]:
    """Extract normalized NCBI Gene identifiers from PubTator output."""

    pairs: list[dict[str, str]] = []

    for document in iter_bioc_documents(payload):
        pmid = str(document.get("id", "")).strip()

        if not pmid:
            continue

        for passage in document.get(
            "passages",
            [],
        ) or []:
            for annotation in passage.get(
                "annotations",
                [],
            ) or []:
                infons = (
                    annotation.get("infons", {})
                    or {}
                )

                concept_type = str(
                    infons.get("type")
                    or infons.get("Type")
                    or ""
                ).lower()

                if (
                    concept_type
                    and concept_type != "gene"
                ):
                    continue

                gene_id = annotation_gene_id(
                    annotation
                )

                if gene_id is None:
                    continue

                pairs.append(
                    {
                        "pmid": pmid,
                        "gene_id": gene_id,
                        "mention": str(
                            annotation.get("text", "")
                        ).strip(),
                    }
                )

    return pairs


def parse_ncbi_gene_info(
    path: Path,
) -> pd.DataFrame:
    """Map human NCBI Gene identifiers to official gene symbols."""

    rows: list[dict[str, str]] = []

    with gzip.open(
        path,
        "rt",
        encoding="utf-8",
    ) as handle:
        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )

        for row in reader:
            if row.get("#tax_id") != "9606":
                continue

            rows.append(
                {
                    "gene_id": row["GeneID"],
                    "official_symbol": (
                        row["Symbol"].upper()
                    ),
                    "description": row.get(
                        "description",
                        "",
                    ),
                    "gene_type_ncbi": row.get(
                        "type_of_gene",
                        "",
                    ),
                }
            )

    return pd.DataFrame(rows)


def normalize_text(text: str) -> str:
    return re.sub(
        r"\s+",
        " ",
        text.lower(),
    ).strip()


def contains_any(
    text: str,
    terms: set[str],
) -> bool:
    return any(term in text for term in terms)


def classify_article_context(
    article: pd.Series,
) -> dict[str, Any]:
    """Assign an evidence class using title and abstract language."""

    title = str(article.get("title", ""))
    abstract = str(article.get("abstract", ""))
    text = normalize_text(f"{title}. {abstract}")

    publication_types = {
        item.strip().lower()
        for item in str(
            article.get(
                "publication_types",
                "",
            )
        ).split("|")
        if item.strip()
    }

    is_review = bool(
        publication_types.intersection(
            REVIEW_PUBLICATION_TYPES
        )
    )

    has_increase = contains_any(
        text,
        INCREASE_TERMS,
    )
    has_decrease = contains_any(
        text,
        DECREASE_TERMS,
    )
    has_benefit = contains_any(
        text,
        BENEFIT_TERMS,
    )
    has_harm = contains_any(
        text,
        HARM_TERMS,
    )
    has_mechanism = contains_any(
        text,
        MECHANISM_TERMS,
    )
    has_model = contains_any(
        text,
        MODEL_TERMS,
    )

    if (
        has_increase
        and has_benefit
        and not is_review
    ):
        evidence_class = "increase_beneficial"

    elif (
        has_decrease
        and has_harm
        and not is_review
    ):
        evidence_class = "decrease_harmful"

    elif (
        has_decrease
        and has_benefit
        and not is_review
    ):
        evidence_class = "decrease_beneficial"

    elif (
        has_increase
        and has_harm
        and not is_review
    ):
        evidence_class = "increase_harmful"

    elif (
        has_mechanism
        and has_model
        and not is_review
    ):
        evidence_class = "mechanistic_primary"

    elif is_review:
        evidence_class = "review"

    else:
        evidence_class = "mention_only"

    return {
        "is_review": is_review,
        "has_model_terms": has_model,
        "has_increase_terms": has_increase,
        "has_decrease_terms": has_decrease,
        "has_benefit_terms": has_benefit,
        "has_harm_terms": has_harm,
        "evidence_class": evidence_class,
    }


def build_gene_article_table(
    pairs: pd.DataFrame,
    articles: pd.DataFrame,
    gene_info: pd.DataFrame,
    downward_symbols: set[str],
) -> pd.DataFrame:
    """Create the gene-publication evidence table."""

    annotated = (
        pairs.merge(
            gene_info,
            on="gene_id",
            how="left",
        )
        .dropna(subset=["official_symbol"])
        .merge(
            articles,
            on="pmid",
            how="left",
        )
    )

    annotated["official_symbol"] = (
        annotated["official_symbol"]
        .astype(str)
        .str.upper()
    )

    annotated = annotated.loc[
        annotated["official_symbol"].isin(
            downward_symbols
        )
    ].copy()

    if annotated.empty:
        return annotated

    context = annotated.apply(
        classify_article_context,
        axis=1,
    )

    context_frame = pd.DataFrame(
        context.tolist(),
        index=annotated.index,
    )

    annotated = pd.concat(
        [
            annotated,
            context_frame,
        ],
        axis=1,
    )

    return annotated.drop_duplicates(
        [
            "official_symbol",
            "pmid",
            "evidence_class",
        ]
    )


def auto_ich_score(
    group: pd.DataFrame,
) -> tuple[int, str, str]:
    """Calculate the preliminary article-derived ICH score."""

    counts = group[
        "evidence_class"
    ].value_counts()

    restorative = int(
        counts.get(
            "increase_beneficial",
            0,
        )
    ) + int(
        counts.get(
            "decrease_harmful",
            0,
        )
    )

    contraindicating = int(
        counts.get(
            "decrease_beneficial",
            0,
        )
    ) + int(
        counts.get(
            "increase_harmful",
            0,
        )
    )

    mechanistic = int(
        counts.get(
            "mechanistic_primary",
            0,
        )
    )

    primary_mentions = int(
        (~group["is_review"]).sum()
    )

    total_mentions = int(
        group["pmid"].nunique()
    )

    if restorative >= 2:
        score = 5
    elif restorative == 1:
        score = 4
    elif (
        contraindicating > 0
        or mechanistic > 0
    ):
        score = 3
    elif primary_mentions > 0:
        score = 2
    elif total_mentions > 0:
        score = 1
    else:
        score = 0

    if (
        restorative > 0
        and contraindicating == 0
    ):
        gate = "restoration_supported"

    elif (
        contraindicating > 0
        and restorative == 0
    ):
        gate = "restoration_contraindicated"

    elif (
        restorative > 0
        and contraindicating > 0
    ):
        gate = "conflicting_direction"

    else:
        gate = "direction_unclear"

    rationale = (
        f"{restorative} restorative-direction primary article(s); "
        f"{contraindicating} contraindicating-direction article(s); "
        f"{mechanistic} mechanistic primary article(s); "
        f"{total_mentions} total annotated article(s)."
    )

    return score, gate, rationale


def aggregate_literature_evidence(
    all_genes: pd.DataFrame,
    gene_articles: pd.DataFrame,
) -> pd.DataFrame:
    """Aggregate the article-level classifications for every gene."""

    if gene_articles.empty:
        result = all_genes.copy()

        count_columns = [
            "ich_article_count",
            "primary_article_count",
            "review_article_count",
            "increase_beneficial_count",
            "decrease_beneficial_count",
            "increase_harmful_count",
            "decrease_harmful_count",
            "mechanistic_primary_count",
            "auto_ich_score",
        ]

        for column in count_columns:
            result[column] = 0

        result["auto_direction_gate"] = (
            "no_ich_literature_detected"
        )
        result["auto_evidence_rationale"] = (
            "No annotated ICH articles detected."
        )
        result["top_pmids"] = ""
        result["top_article_titles"] = ""

        return result

    rows: list[dict[str, Any]] = []

    for gene, group in gene_articles.groupby(
        "official_symbol"
    ):
        score, gate, rationale = auto_ich_score(
            group
        )

        counts = group[
            "evidence_class"
        ].value_counts()

        evidence_order = {
            "increase_beneficial": 0,
            "decrease_harmful": 1,
            "decrease_beneficial": 2,
            "increase_harmful": 3,
            "mechanistic_primary": 4,
            "mention_only": 5,
            "review": 6,
        }

        top_papers = (
            group.assign(
                _evidence_order=(
                    group["evidence_class"]
                    .map(evidence_order)
                    .fillna(9)
                )
            )
            .sort_values(
                [
                    "_evidence_order",
                    "year",
                ],
                ascending=[
                    True,
                    False,
                ],
            )
            .drop_duplicates("pmid")
            .head(5)
        )

        rows.append(
            {
                "Gene": gene,
                "ich_article_count": int(
                    group["pmid"].nunique()
                ),
                "primary_article_count": int(
                    group.loc[
                        ~group["is_review"],
                        "pmid",
                    ].nunique()
                ),
                "review_article_count": int(
                    group.loc[
                        group["is_review"],
                        "pmid",
                    ].nunique()
                ),
                "increase_beneficial_count": int(
                    counts.get(
                        "increase_beneficial",
                        0,
                    )
                ),
                "decrease_beneficial_count": int(
                    counts.get(
                        "decrease_beneficial",
                        0,
                    )
                ),
                "increase_harmful_count": int(
                    counts.get(
                        "increase_harmful",
                        0,
                    )
                ),
                "decrease_harmful_count": int(
                    counts.get(
                        "decrease_harmful",
                        0,
                    )
                ),
                "mechanistic_primary_count": int(
                    counts.get(
                        "mechanistic_primary",
                        0,
                    )
                ),
                "auto_ich_score": score,
                "auto_direction_gate": gate,
                "auto_evidence_rationale": (
                    rationale
                ),
                "top_pmids": "|".join(
                    top_papers[
                        "pmid"
                    ].astype(str)
                ),
                "top_article_titles": " || ".join(
                    top_papers[
                        "title"
                    ].fillna("").astype(str)
                ),
            }
        )

    evidence = pd.DataFrame(rows)

    result = all_genes.merge(
        evidence,
        on="Gene",
        how="left",
    )

    count_columns = [
        "ich_article_count",
        "primary_article_count",
        "review_article_count",
        "increase_beneficial_count",
        "decrease_beneficial_count",
        "increase_harmful_count",
        "decrease_harmful_count",
        "mechanistic_primary_count",
        "auto_ich_score",
    ]

    for column in count_columns:
        result[column] = (
            result[column]
            .fillna(0)
            .astype(int)
        )

    result["auto_direction_gate"] = (
        result["auto_direction_gate"]
        .fillna(
            "no_ich_literature_detected"
        )
    )

    result["auto_evidence_rationale"] = (
        result["auto_evidence_rationale"]
        .fillna(
            "No annotated ICH articles detected."
        )
    )

    result["top_pmids"] = result[
        "top_pmids"
    ].fillna("")

    result["top_article_titles"] = result[
        "top_article_titles"
    ].fillna("")

    return result


def retrieve_articles(
    client: NcbiClient,
    pmids: Sequence[str],
    cache_dir: Path,
) -> pd.DataFrame:
    """Download PubMed abstracts in batches and cache the XML files."""

    records: list[dict[str, Any]] = []

    batches = list(chunks(list(pmids), 200))

    for batch_number, batch in enumerate(
        tqdm(
            batches,
            desc="Retrieving PubMed abstracts",
        )
    ):
        cache_path = (
            cache_dir
            / f"pubmed_{batch_number:05d}.xml"
        )

        if cache_path.exists():
            content = cache_path.read_bytes()
        else:
            content = client.fetch_pubmed_xml(
                batch
            )
            cache_path.parent.mkdir(
                parents=True,
                exist_ok=True,
            )
            cache_path.write_bytes(content)

        records.extend(
            parse_pubmed_articles(content)
        )

    return (
        pd.DataFrame(records)
        .drop_duplicates("pmid")
    )


def retrieve_gene_annotations(
    client: NcbiClient,
    pmids: Sequence[str],
    cache_dir: Path,
) -> pd.DataFrame:
    """Download PubTator gene annotations in batches."""

    pairs: list[dict[str, str]] = []
    batches = list(chunks(list(pmids), 100))

    for batch_number, batch in enumerate(
        tqdm(
            batches,
            desc="Retrieving PubTator annotations",
        )
    ):
        cache_path = (
            cache_dir
            / f"pubtator_{batch_number:05d}.json"
        )

        if cache_path.exists():
            payload = json.loads(
                cache_path.read_text(
                    encoding="utf-8"
                )
            )
        else:
            payload = (
                client.fetch_pubtator_biocjson(
                    batch
                )
            )
            cache_path.parent.mkdir(
                parents=True,
                exist_ok=True,
            )
            cache_path.write_text(
                json.dumps(payload),
                encoding="utf-8",
            )

        pairs.extend(
            parse_pubtator_gene_pairs(payload)
        )

    if not pairs:
        return pd.DataFrame(
            columns=[
                "pmid",
                "gene_id",
                "mention",
            ]
        )

    return (
        pd.DataFrame(pairs)
        .drop_duplicates()
    )


def run_literature_screen(
    downward_genes: pd.DataFrame,
    output_dir: Path,
    email: str,
    api_key: str | None = None,
    query: str = DEFAULT_ICH_QUERY,
    gene_info_url: str = NCBI_GENE_INFO_URL,
) -> pd.DataFrame:
    """Run the full PubMed/PubTator screen for the age-declining genes."""

    output_dir.mkdir(parents=True, exist_ok=True)

    cache_dir = output_dir / "cache"
    raw_dir = output_dir / "raw"

    cache_dir.mkdir(parents=True, exist_ok=True)
    raw_dir.mkdir(parents=True, exist_ok=True)

    settings = ApiSettings(
        email=email,
        api_key=api_key,
    )

    client = NcbiClient(
        settings=settings,
        cache_dir=cache_dir,
    )

    pmid_path = cache_dir / "ich_pmids.json"

    if pmid_path.exists():
        pmids = json.loads(
            pmid_path.read_text(
                encoding="utf-8"
            )
        )
    else:
        pmids = client.search_pubmed(query)
        pmid_path.write_text(
            json.dumps(pmids),
            encoding="utf-8",
        )

    articles_path = raw_dir / "ich_articles.csv"

    if articles_path.exists():
        articles = pd.read_csv(
            articles_path,
            dtype={"pmid": str},
        )
    else:
        articles = retrieve_articles(
            client=client,
            pmids=pmids,
            cache_dir=cache_dir / "pubmed",
        )
        articles.to_csv(
            articles_path,
            index=False,
        )

    annotations_path = (
        raw_dir
        / "pubtator_gene_annotations.csv"
    )

    if annotations_path.exists():
        pairs = pd.read_csv(
            annotations_path,
            dtype={
                "pmid": str,
                "gene_id": str,
            },
        )
    else:
        pairs = retrieve_gene_annotations(
            client=client,
            pmids=pmids,
            cache_dir=cache_dir / "pubtator",
        )
        pairs.to_csv(
            annotations_path,
            index=False,
        )

    gene_info_path = (
        raw_dir
        / "Homo_sapiens.gene_info.gz"
    )

    client.download(
        gene_info_url,
        gene_info_path,
    )

    gene_info = parse_ncbi_gene_info(
        gene_info_path
    )

    gene_articles = build_gene_article_table(
        pairs=pairs,
        articles=articles,
        gene_info=gene_info,
        downward_symbols=set(
            downward_genes["Gene"]
        ),
    )

    gene_articles.to_csv(
        output_dir / "gene_article_evidence.csv",
        index=False,
    )

    screened = aggregate_literature_evidence(
        all_genes=downward_genes,
        gene_articles=gene_articles,
    )

    screened.to_csv(
        output_dir
        / "all_downward_genes_literature_screened.csv",
        index=False,
    )

    return screened
