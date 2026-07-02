# Case 5: Uploaded Result Interpretation

## Input

```text
杩欐槸鎴戞湰鍦拌窇鍑烘潵鐨?marker_gene.csv銆丏EG.csv銆乪nrichment.csv 鍜?cellchat_results.csv锛岃甯垜鍒ゆ柇缁撴灉鏄惁鍚堢悊锛屽苟鐢熸垚璁烘枃 Results 鑽夌銆?```

## Expected Output

- Switch from code generation to result quality review and biological interpretation.
- Check marker support, differential cell groups, key genes, enrichment terms, and CellChat changes.
- Separate observation, statistical inference, mechanism hypothesis, limitation, and validation suggestion.
- Draft cautious Results text without inventing table values.

## Expected Triggered Files

- `agents/05_result_quality_checker.md`
- `agents/06_biological_interpreter.md`
- `agents/07_report_writer.md`
- `agents/08_safety_and_limitation_checker.md`
- `references/output_file_checklist.md`
