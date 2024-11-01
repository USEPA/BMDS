from typing import ClassVar

d_gamma: dich_model
d_hill: dich_model
d_logistic: dich_model
d_loglogistic: dich_model
d_logprobit: dich_model
d_multistage: dich_model
d_probit: dich_model
d_qlinear: dich_model
d_weibull: dich_model
exp_3: cont_model
exp_5: cont_model
funl: cont_model
generic: cont_model
hill: cont_model
log_normal: distribution
nctr: nested_model
nlogistic: nested_model
normal: distribution
normal_ncv: distribution
polynomial: cont_model
power: cont_model

class BMDSMA_results:
    BMD: list[float]
    BMDL: list[float]
    BMDL_MA: float
    BMDU: list[float]
    BMDU_MA: float
    BMD_MA: float
    ebLower: list[float]
    ebUpper: list[float]
    def __init__(self) -> None: ...

class BMDS_results:
    AIC: float
    BIC_equiv: float
    BMD: float
    BMDL: float
    BMDU: float
    bounded: list[bool]
    chisq: float
    lowerConf: list[float]
    slopeFactor: float
    stdErr: list[float]
    upperConf: list[float]
    validResult: bool
    def __init__(self) -> None: ...
    def setSlopeFactor(self, arg0: float) -> None: ...

class cont_model:
    __members__: ClassVar[dict] = ...  # read-only
    __entries: ClassVar[dict] = ...
    exp_3: ClassVar[cont_model] = ...
    exp_5: ClassVar[cont_model] = ...
    funl: ClassVar[cont_model] = ...
    generic: ClassVar[cont_model] = ...
    hill: ClassVar[cont_model] = ...
    polynomial: ClassVar[cont_model] = ...
    power: ClassVar[cont_model] = ...
    def __init__(self, value: int) -> None: ...
    def __and__(self, other: object) -> object: ...
    def __eq__(self, other: object) -> bool: ...
    def __ge__(self, other: object) -> bool: ...
    def __gt__(self, other: object) -> bool: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __int__(self) -> int: ...
    def __invert__(self) -> object: ...
    def __le__(self, other: object) -> bool: ...
    def __lt__(self, other: object) -> bool: ...
    def __ne__(self, other: object) -> bool: ...
    def __or__(self, other: object) -> object: ...
    def __rand__(self, other: object) -> object: ...
    def __ror__(self, other: object) -> object: ...
    def __rxor__(self, other: object) -> object: ...
    def __xor__(self, other: object) -> object: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class continuous_AOD:
    AIC: list[float]
    LL: list[float]
    TOI: testsOfInterest
    addConst: float
    nParms: list[int]
    def __init__(self) -> None: ...

class continuous_GOF:
    calcMean: list[float]
    calcSD: list[float]
    dose: list[float]
    ebLower: list[float]
    ebUpper: list[float]
    estMean: list[float]
    estSD: list[float]
    n: int
    obsMean: list[float]
    obsSD: list[float]
    res: list[float]
    size: list[float]
    def __init__(self) -> None: ...

class dich_model:
    __members__: ClassVar[dict] = ...  # read-only
    __entries: ClassVar[dict] = ...
    d_gamma: ClassVar[dich_model] = ...
    d_hill: ClassVar[dich_model] = ...
    d_logistic: ClassVar[dich_model] = ...
    d_loglogistic: ClassVar[dich_model] = ...
    d_logprobit: ClassVar[dich_model] = ...
    d_multistage: ClassVar[dich_model] = ...
    d_probit: ClassVar[dich_model] = ...
    d_qlinear: ClassVar[dich_model] = ...
    d_weibull: ClassVar[dich_model] = ...
    def __init__(self, value: int) -> None: ...
    def __and__(self, other: object) -> object: ...
    def __eq__(self, other: object) -> bool: ...
    def __ge__(self, other: object) -> bool: ...
    def __gt__(self, other: object) -> bool: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __int__(self) -> int: ...
    def __invert__(self) -> object: ...
    def __le__(self, other: object) -> bool: ...
    def __lt__(self, other: object) -> bool: ...
    def __ne__(self, other: object) -> bool: ...
    def __or__(self, other: object) -> object: ...
    def __rand__(self, other: object) -> object: ...
    def __ror__(self, other: object) -> object: ...
    def __rxor__(self, other: object) -> object: ...
    def __xor__(self, other: object) -> object: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class dicho_AOD:
    devFit: float
    devRed: float
    dfFit: int
    dfRed: int
    fittedLL: float
    fullLL: float
    nFit: int
    nFull: int
    nRed: int
    pvFit: float
    pvRed: float
    redLL: float
    def __init__(self) -> None: ...

class dichotomous_GOF:
    df: float
    ebLower: list[float]
    ebUpper: list[float]
    expected: list[float]
    n: int
    p_value: float
    residual: list[float]
    test_statistic: float
    def __init__(self) -> None: ...

class distribution:
    __members__: ClassVar[dict] = ...  # read-only
    __entries: ClassVar[dict] = ...
    log_normal: ClassVar[distribution] = ...
    normal: ClassVar[distribution] = ...
    normal_ncv: ClassVar[distribution] = ...
    def __init__(self, value: int) -> None: ...
    def __and__(self, other: object) -> object: ...
    def __eq__(self, other: object) -> bool: ...
    def __ge__(self, other: object) -> bool: ...
    def __gt__(self, other: object) -> bool: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __int__(self) -> int: ...
    def __invert__(self) -> object: ...
    def __le__(self, other: object) -> bool: ...
    def __lt__(self, other: object) -> bool: ...
    def __ne__(self, other: object) -> bool: ...
    def __or__(self, other: object) -> object: ...
    def __rand__(self, other: object) -> object: ...
    def __ror__(self, other: object) -> object: ...
    def __rxor__(self, other: object) -> object: ...
    def __xor__(self, other: object) -> object: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class nestedBootstrap:
    pVal: list[float]
    perc50: list[float]
    perc90: list[float]
    perc95: list[float]
    perc99: list[float]
    def __init__(self) -> None: ...

class nestedLitterData:
    LSC: list[float]
    SR: list[float]
    dose: list[float]
    estProb: list[float]
    expected: list[float]
    litterSize: list[float]
    observed: list[int]
    def __init__(self) -> None: ...

class nestedReducedData:
    dose: list[float]
    lowerConf: list[float]
    propAffect: list[float]
    upperConf: list[float]
    def __init__(self) -> None: ...

class nestedSRData:
    avgAbsSR: float
    avgSR: float
    maxAbsSR: float
    maxSR: float
    minAbsSR: float
    minSR: float
    def __init__(self) -> None: ...

class nested_model:
    __members__: ClassVar[dict] = ...  # read-only
    __entries: ClassVar[dict] = ...
    nctr: ClassVar[nested_model] = ...
    nlogistic: ClassVar[nested_model] = ...
    def __init__(self, value: int) -> None: ...
    def __and__(self, other: object) -> object: ...
    def __eq__(self, other: object) -> bool: ...
    def __ge__(self, other: object) -> bool: ...
    def __gt__(self, other: object) -> bool: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __int__(self) -> int: ...
    def __invert__(self) -> object: ...
    def __le__(self, other: object) -> bool: ...
    def __lt__(self, other: object) -> bool: ...
    def __ne__(self, other: object) -> bool: ...
    def __or__(self, other: object) -> object: ...
    def __rand__(self, other: object) -> object: ...
    def __ror__(self, other: object) -> object: ...
    def __rxor__(self, other: object) -> object: ...
    def __xor__(self, other: object) -> object: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class python_continuous_analysis:
    BMD_type: int
    BMR: float
    Y: list[float]
    alpha: float
    burnin: int
    degree: int
    detectAdvDir: bool
    disttype: int
    doses: list[float]
    isIncreasing: bool
    model: cont_model
    n: int
    n_group: list[float]
    parms: int
    prior: list[float]
    prior_cols: int
    restricted: bool
    samples: int
    sd: list[float]
    suff_stat: bool
    tail_prob: float
    transform_dose: int
    def __init__(self) -> None: ...

class python_continuous_model_result:
    aod: continuous_AOD
    bmd: float
    bmd_dist: list[float]
    bmdsRes: BMDS_results
    cov: list[float]
    dist: int
    dist_numE: int
    gof: continuous_GOF
    max: float
    model: int
    model_df: float
    nparms: int
    parms: list[float]
    total_df: float
    def __init__(self) -> None: ...

class python_dichotomousMA_analysis:
    actual_parms: list[int]
    modelPriors: list[float]
    models: list[int]
    nmodels: int
    nparms: list[int]
    prior_cols: list[int]
    priors: list[list[float]]
    pyDA: python_dichotomous_analysis
    def __init__(self) -> None: ...

class python_dichotomousMA_result:
    bmd_dist: list[float]
    bmdsRes: BMDSMA_results
    dist_numE: int
    models: list[python_dichotomous_model_result]
    nmodels: int
    post_probs: list[float]
    def __init__(self) -> None: ...

class python_dichotomous_analysis:
    BMD_type: int
    BMR: float
    Y: list[float]
    alpha: float
    burnin: int
    degree: int
    doses: list[float]
    model: int
    n: int
    n_group: list[float]
    parms: int
    prior: list[float]
    prior_cols: int
    samples: int
    def __init__(self) -> None: ...

class python_dichotomous_model_result:
    aod: dicho_AOD
    bmd: float
    bmd_dist: list[float]
    bmdsRes: BMDS_results
    cov: list[float]
    dist_numE: int
    gof: dichotomous_GOF
    gof_chi_sqr_statistic: float
    gof_p_value: float
    max: float
    model: int
    model_df: float
    nparms: int
    parms: list[float]
    total_df: float
    def __init__(self) -> None: ...

class python_multitumor_analysis:
    BMD_type: int
    BMR: float
    alpha: float
    degree: list[int]
    models: list[list[python_dichotomous_analysis]]
    n: list[int]
    ndatasets: int
    nmodels: list[int]
    prior_cols: int
    def __init__(self) -> None: ...

class python_multitumor_result:
    BMD: float
    BMDL: float
    BMDU: float
    combined_LL: float
    combined_LL_const: float
    models: list[list[python_dichotomous_model_result]]
    ndatasets: int
    nmodels: list[int]
    selectedModelIndex: list[int]
    slopeFactor: float
    validResult: bool
    def __init__(self) -> None: ...
    def setSlopeFactor(self, arg0: float) -> None: ...

class python_nested_analysis:
    BMD_type: int
    BMR: float
    ILC_type: int
    LSC_type: int
    alpha: float
    doses: list[float]
    estBackground: bool
    incidence: list[float]
    iterations: int
    litterSize: list[float]
    lsc: list[float]
    model: nested_model
    numBootRuns: int
    parms: int
    prior: list[float]
    prior_cols: int
    seed: int
    def __init__(self) -> None: ...

class python_nested_result:
    LL: float
    bmd: float
    bmdsRes: BMDS_results
    boot: nestedBootstrap
    combPVal: float
    cov: list[float]
    fixedLSC: float
    litter: nestedLitterData
    model: nested_model
    model_df: float
    nparms: int
    parms: list[float]
    reduced: nestedReducedData
    srData: nestedSRData
    validResult: bool
    def __init__(self) -> None: ...

class test_struct:
    BMD: float
    doses: list[float]
    n: int
    validResult: bool
    def __init__(self) -> None: ...

class testsOfInterest:
    DF: list[float]
    llRatio: list[float]
    pVal: list[float]
    def __init__(self) -> None: ...

def pythonBMDSCont(
    python_continuous_analysis: python_continuous_analysis,
    python_continuous_model_result: python_continuous_model_result,
) -> None: ...
def pythonBMDSDicho(
    python_dichotomous_analysis: python_dichotomous_analysis,
    python_dichotomous_model_result: python_dichotomous_model_result,
) -> None: ...
def pythonBMDSDichoMA(
    python_dichotomousMA_analysis: python_dichotomousMA_analysis,
    python_dichotomousMA_result: python_dichotomousMA_result,
) -> None: ...
def pythonBMDSMultitumor(
    python_multitumor_analysis: python_multitumor_analysis,
    python_multitumor_result: python_multitumor_result,
) -> None: ...
def pythonBMDSNested(
    python_nested_analysis: python_nested_analysis, python_nested_result: python_nested_result
) -> None: ...
def version() -> str: ...
