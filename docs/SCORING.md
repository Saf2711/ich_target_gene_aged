# Scoring definitions

## Age score

| Discovery tier | Age score |
|---|---:|
| High-confidence age decline | 5 |
| Moderate age-decline evidence | 4 |
| Directional age decline | 3 |
| Not prioritized | 0 |

## Automated ICH evidence score

| Score | Evidence rule |
|---:|---|
| 5 | At least two restorative-direction primary articles |
| 4 | One restorative-direction primary article |
| 3 | Mechanistic primary evidence or contraindicating-direction intervention |
| 2 | Primary article mention |
| 1 | Review-only evidence |
| 0 | No annotated ICH article |

Restorative-direction evidence includes:

```text
increase + benefit
decrease + harm
```

Contraindicating-direction evidence includes:

```text
decrease + benefit
increase + harm
```

## Payload score

| Coding sequence | Score |
|---:|---:|
| ≤800 bp | 5 |
| 801–1,500 bp | 4 |
| 1,501–2,250 bp | 3 |
| 2,251–3,000 bp | 2 |
| 3,001–4,000 bp | 1 |
| >4,000 bp or unavailable | 0 |

## Final score

Default weights:

```text
Age evidence: 0.20
ICH evidence: 0.60
Payload suitability: 0.20
```

Formula:

```text
20 × [(age × 0.20) + (ICH × 0.60) + (payload × 0.20)]
```

The maximum is 100.
