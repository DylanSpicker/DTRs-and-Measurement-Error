# Dynamic Treatment Regimes + Measurement Error
The code contained within this repository has been created for the facilitation of my Master's Thesis. The results, in conjunction with theoretical development, are currently being written-up for the Masters project. Except where referenced in the code, the work is entirely mine, and it should be considered a WIP.

## Theoretical Background
My research focuses on fitting Dynamic Treatment Regimes (DTRs) in the pressence of measurement error. The following is a quick overview on the content, so as to familiarize with my terminology, etc. but for a comprehensive overview of the subject areas, consider the references section.

### DTRs
A DTR, at its core, is a manner of modelling a decision making process rigorously. Typically, this process is assumed to be the treatment process for a disease. The idea is to encode an iterative process in a statistical model which can be readily fit, but which can accomodate diverse structures. If we assume that there are some set of covariates (denoted X) that we care about (i.e. blood pressure, age, weight, etc.), and we are focused on affecting some health outcome (denoted Y), then a DTR revolves arrounding finding the optimal treatment decision (denoted A), to maximize the outcome Y, given all the covariates in X. In my work, I tend to assume that A is binary (either 0 or 1), which may correspond to a control treatment and an experimental treatment, for instance. 

Overall then, a DTR asks "When a patient presents with covariate measures X, will Y be higher if A=1 or if A=0?" In order to model this we consider defining a treatment model, a treatment-free model, and a blip model. The treatment-free model captures the effect that the covariates have on "Y" in the absence of treatment; the blip model captures the added benefit (the blip) that occurs in "Y" from treatment; the treatment model captures the probability that a random individual in our sample was assigned to treatment, given their covariates. So for instance, we may assume that:

```Y = 1 + X + A*(2 - 3*X)```

and suggest that: 

```P(A=1|X=x) = x```

In this case, we would say that the "treatment-free" model is "1+X", that the "blip model" is "2-3*X", and that the treatment model is "x". 

This framework is particularly powerful as it can be extended to multiple decision points. That is, a patient presents with some set of covariates, and a treatment decision is made. This treatment may impact the covariates, and may impact the level Y. After sometime, the patient is seen again, and a second treatment decision can be made. We can do this iteratively for as long as any treatment procedure occurs, capturing a wide array of medical proceedings (and non-medical related decisions). While the estimation of the optimal DTRs is outside of the scope of this brief overview, note that it can be done consistently in the multistage setting, leveraging simple weighted regression techniques.

### Measurement Error
Measurement error captures the idea that often we are unable to observe the true value of some covariate. Instead, we may measure a value that is prone-to-error, which we use in place of the true value. This may come as a consequence of the measurement we use (i.e. blood pressure can only be accurately measured to some tolerance level), through a systematic mechanism (i.e. patients will consistently under report their caloric intake), or due to a proxy of the truth (i.e. measuring the dilation of blood vessels near the surface of the skin to approximate the dilation of internal vessels), but in any event it ensures that standard analysis will tend to be unreliable. 

In general, in my research, I will assume a classical additive measurement error. That is, instead of observing the true covariates X, we observe W = X + U, where U follows some normal distribution with 0 mean and fixed variance. U is assumed to be independent of X. The basic premise of the research is to investigate how the presence of error impacts DTRs, and to attempt to see whether there are easily implementable strategies to overcome any issues which arise.

## Regression Calibration
This entire repository relies on my implementation of Regression Calibration, available at [DylanSpicker/rCalibration](https://github.com/DylanSpicker/rCalibration). Full details are available at the linked repository.

## Repository Structure
The code presented within produces simulations for a number of results which are needed to be shown for the paper. Broadly speaking, the code has been written to preserve as much of the data as possible, separating any relevant analysis from the simulations themselves. I have also avoided using any R packages, outside of the standard library, for no reason other than ease of transfer to other researchers. At present there are five overarching sets of simulations that are included, each outlined below. The data from the process should be accessible in the `.RData` files contained in the relevant directories.

### Basic Error-Prone Simulations
* Directory: `basic-sim-41/`

These simulations simply look at a basic DTR, with measurement error, and seek to demonstrate the issues that arise using a naive analysis when error is introduced. Further, it demonstrates that a simple [regression calibration](https://github.com/DylanSpicker/rCalibration) correction can be applied to account for the bias in the naive estimators. In particular, this is showcased in the most simplistic model where Y = 1 + X + A(1 + X), the treatment model is a simple logistic regression, and everything is correctly specified. 

### Double Robustness
* Directory: `double-robustness-42/`

These simulations attempt to look at the double-robustness of the dWOLS estimation procedure. In particular the simulations show that, even if one of the treatment or treatment-free models is misspecified, so long as the blip is correctly specified, that dWOLS can lead to accurate measurements using regression calibration corrected covariates. Here we leverage 20% replication of error-prone covariates, and simulate Y = 1 - X + exp(X) + A(3 - 2X), with a treatment model specified as H(exp(W1)), where H() denotes the standard expit function. 

### Independence Violations
* Directory: `independence-violations-43/`

These simulations attempt to push the limits of the interaction between logistic regression and regression calibration; we rely on accurate estimates for the propensity scores to generate independence between the covariates and the treatments in the models, and so these simulations investigate the level to which our estimates remain accurate when the independence assumption is violated on account of poorly modeled propensity scores. These simulations vary the treatment model dramatically, over 36 pairs of coefficients for the logistic regression, fixing the marginal probability of treatment to be one of {0.1, 0.5, 0.9}, while altering the effect-size in the treatment-model.

### Predictability of Bias
* Directory: `predictability-of-bias-44/`

These simulations seek to demonstrate that, while the previous naive analysis appears to always attenuate in line with a standard reliability ratio, this will not always be the case in DTRs with ME. In particular, it is possible to have attenuation or amplification of the bias in a naive analysis. At the present, I am still determining which scenarios to include in the write-up to best demonstrate these possibilities; as such, only two of the four scenarios are currently being run (and consequently, only 2 of the scenarios are currently saved in the `.RData`). The other scenarios should run perfectly if uncommented in the code, but it remains to be seen as to whether or not they show anything of value outside of the two scenarios presented.

### Future Treatment
* Directory `future-treatments-45/`

These simulations differ from the others in the repository in that a complete DTR is not fit. Instead, we seek to determine what will happen when assigning future treatments based on estimated blip parameters. In particular, if we assume that we have correctly estimated the blip parameters from a study, and are now looking to use these to treat incoming patients, what happens if the practitioners measurements are still prone to error? Should we attempt to use an adjusted measurement (i.e. regression calibration on the incoming patients), or can we simply use the naive estimate? Further, if we assume that the quantity of replicates available for use is limited (i.e. can only replicate 10% or 50% of all measures), then we can attempt to discern how these replicates should be allocated (i.e. should we randomly select them, or apply them to the observations closest to the decision boundary) to be facilitate optimal treatment.

## References
* Chakraborty, B. & Moodie, E. E. (2013). Statistical methods for dynamic treatment regimes. Springer New York. doi:10.1007/978-1-4614-7428-9
* Wallace, M. P. & Moodie, E. E. M. [Erica E. M.]. (2014, April). Personalizing medicine: a review of adaptive treatment strategies. Pharmacoepidemiology and Drug Safety, 23 (6), 580–585. doi:10.1002/pds.3606
* Wallace, M. P. & Moodie, E. E. M. [Erica E M]. (2015). Doubly-robust dynamic treatment regimen estimation via weighted least squares. Biometrics, 71 (3), 636–644. doi:10.1111/biom.12306
