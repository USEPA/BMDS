# BMDS

<table>
<tbody>
<tr>
    <td>Package</td><td>

[![version](https://img.shields.io/pypi/v/pybmds.svg?label=pybmds%20version&maxAge=3600)](https://pypi.org/project/pybmds/)
[![version](https://img.shields.io/pypi/v/bmds-ui.svg?label=bmds-ui%20version&maxAge=3600)](https://pypi.org/project/bmds-ui/)
    </td>
</tr>
<tr>
    <td>Documentation</td><td>

[![Docs Badge](https://img.shields.io/badge/Latest-online-brightgreen)](https://usepa.github.io/BMDS) [![Read The Docs](https://img.shields.io/badge/Versioned-online-brightgreen)](https://pybmds.readthedocs.io/)
    </td>
</tr>
<tr>
    <td>Website</td><td>

[![Site Badge](https://img.shields.io/badge/BMDS-online-purple)](https://epa.gov/bmds) [![Site Badge](https://img.shields.io/badge/BMDS&nbsp;Online-online-purple)](https://bmdsonline.epa.gov)
    </td>
</tr>

</tbody>
</table>


EPA's Benchmark Dose Software (BMDS) collects and provides easy access to numerous mathematical models that help risk assessors estimate the quantitative relationship between a chemical dose and the test subject’s response.  A specific focus of BMDS is estimating a statistical benchmark dose (BMD). The BMD is a chemical dose or concentration that produces a predetermined change in the response rate of an adverse effect, such as weight loss or tumor incidence. The BMD is a range, rather than a fixed number. For example, the benchmark dose (lower confidence limit) (BMDL) can be regarded as a dose where the observable physical effect is less than the predetermined benchmark response (BMR).

![An example dose response plot an and curve fit](https://github.com/USEPA/BMDS/raw/e89f79388dc3021604e1230ac75e721c12c6bf61/tests/test_pybmds/data/mpl/test_dichotomous_plot.png)

Additional information, documentation, and technical guidance for BMDS is available at [https://www.epa.gov/bmds](https://www.epa.gov/bmds).

This repository contains a low-level C++ library, `bmdscore`, and the `pybmds` Python package for interfacing with `bmdscore` with higher level utilities such as plotting and reporting.

## Credits

The authors would like to thank Dr. Matt Wheeler of NIH/NIEHS for his contributions to many of the `bmdscore` algorithms, and his continued collaboration with the [ToxicR](https://github.com/NIEHS/ToxicR) software.

## Disclaimer

The United States Environmental Protection Agency (EPA) GitHub project code is provided on an "as is" basis and the user assumes responsibility for its use.  EPA has relinquished control of the information and no longer has responsibility to protect the integrity, confidentiality, or availability of the information.  Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or favoring by EPA.  The EPA seal and logo shall not be used in any manner to imply endorsement of any commercial product or activity by EPA or the United States Government.
