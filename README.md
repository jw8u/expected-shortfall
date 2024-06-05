# ES

Goal is to replicate some results from a paper via simulations. In particular, the $\hat{\text{ES}}$ (adjusted tail-based normal approximation ES estimator) and the AA (Arithmetic Average ES estimator) are compared. 

## Usage
```
Rscript ES.R
```

## Output
*print_t(M)* looks at the *t*-distribution with *df = 3.5, 5, 8*, *n = 250, 500* and $\beta = 0.99, 0.995$  and prints out the MSE, variance, and bias of both the $\hat{\text{ES}}$ and AA ES estimators through a Monte Carlo simulation of *M* replications.
*print_gpd(M)*, *print_weibull(M)*, *print_t2(M)*, and *print_gpd2(M)* are defined similarly for other parameters/distributions.

The output for the above functions for *M = 2500* can be re-organized into the tables in the following section. Please note that all numerical results will be identical upon each running of the script - this is because seed is set to 001 before each replication.

For full explanations of methodology and results, please see *Report.pdf*

## Numerical Results

**t-distribution, n = 250**
|        |                   | MSE $\beta = 0.99$ | MSE $\beta = 0.995$ | MSE $\beta = 0.99$ | MSE $\beta = 0.995$ | MSE $\beta = 0.99$ | MSE $\beta = 0.995$ |
|--------|-------------------|--------------------|---------------------|--------------------|---------------------|--------------------|---------------------|
| df=3.5 | $\hat{\text{ES}}$ | 3.513701           | 6.918576            | 3.349879           | 5.931118            | -0.4064011         | -0.9949018          |
|        | AA                | 3.045859           | 6.618685            | 2.757139           | 4.895043            | -0.5383521         | -1.313621           |
| df=5   | $\hat{\text{ES}}$ | 0.9119978          | 1.772814            | 0.8458146          | 1.465247            | -0.2579177         | -0.5551157          |
|        | AA                | 0.896156           | 1.88101             | 0.7838646          | 1.292381            | -0.3355666         | -0.7675588          |
| df=8   | $\hat{\text{ES}}$ | 0.3733118          | 0.6882679           | 0.3424361          | 0.5728771           | -0.1761041         | -0.3400294          |
|        | AA                | 0.3821088          | 0.7564698           | 0.3286057          | 0.5127364           | -0.2315912         | -0.4939012          |

**t-distribution, n = 500**
|        |                   | MSE $\beta = 0.99$ | MSE $\beta = 0.995$ | MSE $\beta = 0.99$ | MSE $\beta = 0.995$ | MSE $\beta = 0.99$ | MSE $\beta = 0.995$ |
|--------|-------------------|--------------------|---------------------|--------------------|---------------------|--------------------|---------------------|
| df=3.5 | $\hat{\text{ES}}$ | 2.32565            | 4.622969            | 2.311198           | 4.324478            | 0.1239996          | 0.547924            |
|        | AA                | 1.683792           | 4.569174            | 1.514917           | 4.20693             | 0.4116804          | 0.6032632           |
| df=5   | $\hat{\text{ES}}$ | 0.548986           | 1.095874            | 0.5352389          | 0.9876675           | 0.1181574          | 0.329547            |
|        | AA                | 0.4980796          | 1.152373            | 0.42303            | 1.005189            | 0.2742606          | 0.3841695           |
| df=8   | $\hat{\text{ES}}$ | 0.1927528          | 0.3692167           | 0.1858955          | 0.3318638           | -0.08325665        | -0.1936121          |
|        | AA                | 0.1933235          | 0.4038955           | 0.1597503          | 0.3446732           | -0.183404          | -0.2436394          |

**weibull(shape, scale=1)-distribution, n = 250**
|           |                   | MSE $\beta = 0.99$ | MSE $\beta = 0.995$ | MSE $\beta = 0.99$ | MSE $\beta = 0.995$ | MSE $\beta = 0.99$ | MSE $\beta = 0.995$ |
|-----------|-------------------|--------------------|---------------------|--------------------|---------------------|--------------------|---------------------|
| shape=0.6 | $\hat{\text{ES}}$ | 21.13662           | 40.22232            | 19.36554           | 32.91464            | -1.333727          | -2.705706           |
|           | AA                | 20.8001            | 42.5514             | 18.02868           | 28.79762            | -1.666921          | -3.710162           |
| shape=0.9 | $\hat{\text{ES}}$ | 1.249634           | 2.257493            | 1.12784            | 1.846376            | -0.3496358         | -0.6417602          |
|           | AA                | 1.306666           | 2.532453            | 1.100847           | 1.638437            | -0.4541577         | -0.9458706          |
| shape=1.4 | $\hat{\text{ES}}$ | 0.127255           | 0.2189877           | 0.1145921          | 0.1810451           | -0.1127329         | -0.1949746          |
|           | AA                | 0.1376682          | 0.2554103           | 0.1137437          | 0.1608685           | -0.1548224         | -0.3075811          |

**weibull(shape, scale=1)-distribution, n = 500**
|           |                   | MSE $\beta = 0.99$ | MSE $\beta = 0.995$ | MSE $\beta = 0.99$ | MSE $\beta = 0.995$ | MSE $\beta = 0.99$ | MSE $\beta = 0.995$ |
|-----------|-------------------|--------------------|---------------------|--------------------|---------------------|--------------------|---------------------|
| shape=0.6 | $\hat{\text{ES}}$ | 11.74906           | 23.06104            | 11.26571           | 20.20684            | -0.6984654         | -1.691829           |
|           | AA                | 11.25699           | 24.61849            | 9.428502           | 21.14317            | -1.353611          | -1.866489           |
| shape=0.9 | $\hat{\text{ES}}$ | 0.6679525          | 1.240388            | 0.635357           | 1.09569             | -0.1812447         | -0.3809666          |
|           | AA                | 0.7109973          | 1.378274            | 0.580438           | 1.155482            | -0.361651          | -0.4724979          |
| shape=1.4 | $\hat{\text{ES}}$ | 0.06647344         | 0.1170228           | 0.06338532         | 0.1050339           | -0.05579845        | -0.1096854          |
|           | AA                | 0.0748008          | 0.1327834           | 0.06002494         | 0.1094981           | -0.1216547         | -0.1527384          |

**Generalized Pareto Distribution ($\xi = 0.3,0.2,0.1$), n = 250**
|             |                   | MSE $\beta = 0.99$ | MSE $\beta = 0.995$ | MSE $\beta = 0.99$ | MSE $\beta = 0.995$ | MSE $\beta = 0.99$ | MSE $\beta = 0.995$ |
|-------------|-------------------|--------------------|---------------------|--------------------|---------------------|--------------------|---------------------|
| $\xi = 0.3$ | $\hat{\text{ES}}$ | 38.5636            | 75.88232            | 37.04943           | 65.80364            | -1.236523          | -3.178836           |
|             | AA                | 31.10906           | 69.26451            | 28.36118           | 52.39458            | -1.661093          | -4.109853           |
| $\xi = 0.2$ | $\hat{\text{ES}}$ | 8.775228           | 17.12799            | 8.216385           | 14.30545            | -0.7497529         | -1.681746           |
|             | AA                | 8.017563           | 17.15953            | 7.086977           | 12.14669            | -0.9661373         | -2.240021           |
| $\xi = 0.1$ | $\hat{\text{ES}}$ | 2.298071           | 4.326176            | 2.10109            | 3.535198            | -0.4447709         | -0.890164           |
|             | AA                | 2.290295           | 4.661933            | 1.964367           | 3.119049            | -0.5715881         | -1.242631           |

**Generalized Pareto Distribution ($\xi = 0.3,0.2,0.1$), n = 500**
|             |                   | MSE $\beta = 0.99$ | MSE $\beta = 0.995$ | MSE $\beta = 0.99$ | MSE $\beta = 0.995$ | MSE $\beta = 0.99$ | MSE $\beta = 0.995$ |
|-------------|-------------------|--------------------|---------------------|--------------------|---------------------|--------------------|---------------------|
| $\xi = 0.3$ | $\hat{\text{ES}}$ | 29.75435           | 59.14329            | 29.56326           | 55.50351            | -0.4504673         | -1.913632           |
|             | AA                | 18.02656           | 51.4315             | 16.05085           | 47.03413            | -1.407881          | -2.101471           |
| $\xi = 0.2$ | $\hat{\text{ES}}$ | 5.718911           | 11.47192            | 5.594679           | 10.39187            | -0.3556264         | -1.041253           |
|             | AA                | 4.551092           | 11.32327            | 3.893655           | 9.993395            | -0.8117846         | -1.154933           |
| $\xi = 0.1$ | $\hat{\text{ES}}$ | 1.344728           | 2.624374            | 1.291502           | 2.324945            | -0.2318252         | -0.5480504          |
|             | AA                | 1.288242           | 2.796771            | 1.063535           | 2.38474             | -0.4744814         | -0.6426391          |

## Authors

- [@jw8u](https://www.github.com/jw8u)
