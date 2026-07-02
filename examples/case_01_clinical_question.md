# Case 1: Clinical Research Question

## Input

```text
鎴戞兂鐮旂┒鑶濋鍏宠妭鐐庢粦鑶滅粍缁囦腑宸ㄥ櫖缁嗚優鍜屾垚绾ょ淮缁嗚優鐨勫紓甯稿彉鍖栵紝鎯崇敤 GEO 鍗曠粏鑳炴暟鎹仛鏈哄埗鍒嗘瀽銆?```

## Expected Output

- Identify disease: knee osteoarthritis.
- Identify tissue: synovium.
- Identify key cells: macrophages, fibroblasts, T cells, endothelial cells.
- Recommend modules: dataset search/diagnosis, QC, annotation, cell proportion, DEG, enrichment, CellChat, and pseudotime when justified.
- Ask for GEO accession, organism, comparison design, and sample metadata if missing.

## Expected Triggered Files

- `agents/01_clinical_question_parser.md`
- `agents/03_analysis_plan_generator.md`
- `references/clinical_translation_rules.md`
