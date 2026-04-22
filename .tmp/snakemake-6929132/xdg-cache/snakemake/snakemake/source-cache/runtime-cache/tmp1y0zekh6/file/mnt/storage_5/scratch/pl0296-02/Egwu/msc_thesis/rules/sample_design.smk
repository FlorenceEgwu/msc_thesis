rule sample_design:
    output:
        sample_design_target()
    run:
        workflowSampleDesign.write_sample_design(output[0], SAMPLE_ROWS, SAMPLE_IDS)
