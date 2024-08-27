from pybmds.reporting.styling import Report, write_dataset_table


def test_write_dataset_table(cdataset, cidataset, ddataset, nd_dataset):
    report = Report.build_default()
    write_dataset_table(report, cdataset, True)
    write_dataset_table(report, cdataset, False)
    write_dataset_table(report, cidataset, True)
    write_dataset_table(report, cidataset, False)
    write_dataset_table(report, ddataset, True)
    write_dataset_table(report, ddataset, False)
    write_dataset_table(report, nd_dataset, False)
