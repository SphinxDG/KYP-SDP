# KYP-SDP
This repository contains a structure exploiting algorithm for Robust Control problems.

Many Robust Control problems can be posed as linear matrix equation (LMI) optimization problems of the form

$$
\begin{align}
	&\min_{\lambda \in \mathbb{R}^p, P \in \mathbb{S}^n} ~~ c^\top \lambda - trace (\Sigma P)\\
	&\mathrm{s.t.} ~~
	\begin{pmatrix}
		A & B\\
		I & 0
	\end{pmatrix}^\top 
	\begin{pmatrix}
		0 & P\\
		P & 0
	\end{pmatrix}
	\begin{pmatrix}
		A & B\\
		I & 0
	\end{pmatrix}
	+
	\begin{pmatrix}
		Q(\lambda) & S(\lambda)\\
		S^\top (\lambda) & R(\lambda)
	\end{pmatrix}
	\prec 0, \\
	&\hspace{8mm} N(\lambda) \succ 0,\\
    &\hspace{8mm} P \succ 0.
\end{align}
$$

In this optimization problem, the matrix pair $(A,B)$ defines a linear dynamical system

$$
\dot{x} = A x + Bu,
$$

and is assumed to be controllable. Furthermore, $Q(\cdot)$, $S(\cdot)$, $R(\cdot)$ and $N(\cdot)$ are affine
functions of the optimization variable $\lambda$, i.e., $H(\lambda) := H_0 + \sum_{i = 1}^p \lambda_i H_i$ 
for $H \in \{N, Q, S, R\}$. The optimization problem above has been called KYP-SDP in some literature sources,
because the linear matrix inequality is related to the celebrated KYP-Lemma from Robust Control.

The solver provided in this repository exploits the structure of the optimization problem by eliminating the
matrix variable $P$ with the solution of a Riccati equation. This should be particularly efficient, when the
state dimension $n$ is larger than the number of multipliers $p$.

Our algorithm expects matrices $Q_i \in \mathbb{S}^n$, $S_i \in \mathbb{R}^{n\times m}$, $R_i \in \mathbb{S}^m$, $N_i \in \mathbb{S}^r$,
$A \in \mathbb{R}^{n\times n}$, $B \in \mathbb{R}^{n\times m}$, $\Sigma \in \mathbb{S}^n$ and $c \in \mathbb{R}^p$ for $i = 0,\ldots, p$
as input. Here, $\Sigma$ must be a positive definite matrix and $n,m,p,r \in \mathbb{N}$ must be positive integers.


The algorithm is contained in the file
 - customKYPLMIsolver.m

In addition to the code for the solver, this repository contains a benchmark example for the solver, which is contained
in the files
 - compLibTestCustomSolver.m
 - robustSynthesisCustomSolver.m
 - robustSynthesisYalmip.m

This experiment is based on the library of control systems provided in http://www.complib.de . Here, "compLibTestCustomSolver.m"
is the main script, whereas "robustSynthesisCustomSolver.m" and "robustSynthesisYalmip.m" contain subroutines for calling
different solvers for robust controller synthesis. Furthermore, this experiment comes with the files
 - loadMatricesFromCompLib.m
 - generateTables.m

to extract data from Complib and to generate the table with computation times provided below.



# Simulation example

To demonstrate the structure exploitation of our algorithm, a benchmark of robust controller
synthesis was carried out. To this end, a dynamical system of the structure

$$
    \dot{x} = \mathcal{A} x + \mathcal{B} \begin{pmatrix}
    1 + \delta_1 & & \\
    & \ddots &\\
    & & 1 + \delta_m
    \end{pmatrix}
    u
$$

is considered, where $\delta_1,\ldots,\delta_m\in [-0.25,0.25]$ are uncertain parameters that model 
an uncertainty in the actuators of the system. Our aim is the design of a robust LQR controller, i.e.,
we want to find a control Lyapunov function $V: \mathbb{R}^n \to \mathbb{R}_{\geq 0}, x \mapsto x^\top P^{-1} x$,
such that there exists a controller $u = Kx$, which guarantees the performance 

$$
    V(x) \geq \int_0^\infty x^\top \mathcal{Q} x + u^\top \mathcal{R} u \mathrm{d}x .
$$


By interpreting the signal

$$
    \begin{pmatrix}
    1 + \delta_1 & & \\
    & \ddots &\\
    & & 1 + \delta_m
    \end{pmatrix}
    u
    =
    u + w
$$

as the sum of the input and a disturbance signal $w$ with $w_i := \delta_i u_i$, we can denote the
control system as

$$
    \dot{x} = \mathcal{A} x + \mathcal{B} u + \mathcal{B} w.
$$

Note that $w$ satisfies the quadratic constraint

$$
	\begin{pmatrix}
		u\\
		w
	\end{pmatrix}^\top M(\lambda)^{-1}
	\begin{pmatrix}
		u\\
		w
	\end{pmatrix} \geq 0 
$$

for

$$
    M(\lambda) = diag (\gamma^2\lambda_1,\ldots,\gamma^2\lambda_p,-\lambda_1,\ldots,-\lambda_p)
$$

and all $\lambda \in \mathbb{R}^p_{> 0}$ and $\gamma = 0.25$. Consequently, $V: \mathbb{R}^n \to \mathbb{R}_{\geq 0}, x \mapsto x^\top P^{-1} x$ 
is such a robust control Lyapunov function if it satisfies

$$
	\min_{u\in \mathbb{R}^m}\nabla V(x)^\top (\mathcal{A} x + \mathcal{B} u + \mathcal{B} w) + x^\top \mathcal{Q} x + u^\top \mathcal{R} u \leq 0
$$

for all $x \in \mathbb{R}^n$ and all $w$ for which the quadratic constraint holds. Using (standard) multiplier relaxation and
elimination techniques from Robust Control, we can formulate the linear matrix inequality optimization problem

$$
    \begin{align}
    \min_{\lambda \in \mathbb{R}^p, P \in \mathbb{S}^n} ~&~ c^\top \lambda - trace ( P)\\
    \text{s.t.   } &
    \begin{pmatrix}
		\mathcal{A}^\top & I & C^\top\\
		I & 0 & 0
	\end{pmatrix}^\top
	\begin{pmatrix}
		0 & P\\
		P & 0
	\end{pmatrix}
	\begin{pmatrix}
		\mathcal{A}^\top & I & C^\top\\
		I & 0 & 0
	\end{pmatrix}
	- 
	\begin{pmatrix}
		0 & -I & 0\\
		\mathcal{B}^\top & 0 & I\\
		0 & 0 & -I\\
		\mathcal{B}^\top & 0 & 0
	\end{pmatrix}^\top 
	\begin{pmatrix}
		\mathcal{Q}^{-1} & 0 & &\\
		0 & \mathcal{R}^{-1} & &\\
		& & M_{11}(\lambda) & M_{12}(\lambda)\\
		& & M_{21}(\lambda) & M_{22}(\lambda)
	\end{pmatrix}
	\begin{pmatrix}
		0 & -I & 0\\
		\mathcal{B}^\top & 0 & I\\
		0 & 0 & -I\\
		\mathcal{B}^\top & 0 & 0
	\end{pmatrix} \prec 0\\
	~&~
	\lambda_i > 0 ~~~~~~~ \forall i = 1,\ldots,p\\
	~&~P \succ 0
    \end{align}
$$

for $P$ and $\lambda$. This is a KYP-SDP for robust state feedback controller synthesis.

To benchmark our solver, we extract matrices $(\mathcal{A},\mathcal{B})$ from CompLib and solve the KYP-SDP above for various dynamical systems.

The following table compares computation times (in seconds) of the structure exploiting solver and the three off-the-shelve
solvers Mosek, SeDuMi and LMILab. Inf means that the solver did not terminate after 10^4 seconds. (We appologies that these values do not match
exactly those in our publication. They have been generated in an independent simulation run.)


    Problem      customSolver     Mosek       SeDuMi      LMILab 
    __________    ____________    ________    ________    ________

    {'AC1'   }       0.10389        1.0029      1.5945      1.7903
    {'AC2'   }      0.044277       0.22242     0.25639     0.29568
    {'AC3'   }      0.022823       0.17959     0.19912     0.19895
    {'AC4'   }      0.028846       0.11923     0.14967     0.13055
    {'AC5'   }      0.059092       0.49323     0.54303     0.51754
    {'AC6'   }      0.021965      0.084626      0.1174     0.28713
    {'AC7'   }      0.025528      0.083923      0.1187        1.12
    {'AC8'   }      0.014432       0.07211     0.11523      1.1686
    {'AC9'   }      0.064812      0.076574      0.1373      3.2895
    {'AC10'  }       0.32521        3.7069      20.363         Inf
    {'AC11'  }      0.015928      0.066153    0.091478    0.095954
    {'AC12'  }       0.04962      0.065037    0.095246    0.078521
    {'AC13'  }       0.26896       0.25955      1.1019      7712.6
    {'AC14'  }       0.46766       0.70465      4.6848         Inf
    {'AC15'  }      0.016287      0.063627    0.090293     0.16018
    {'AC16'  }      0.016645      0.065532    0.081684     0.07582
    {'AC17'  }      0.010521      0.065836     0.08138    0.070606
    {'AC18'  }      0.036861      0.074354     0.14215      3.4057
    {'HE1'   }      0.019001      0.065402    0.082946    0.070229
    {'HE2'   }      0.016085      0.065446    0.084281    0.071446
    {'HE3'   }      0.037622      0.068607     0.11243      0.6285
    {'HE4'   }      0.064238      0.068843     0.11381     0.71555
    {'HE5'   }      0.063774      0.069744     0.10878     0.71924
    {'HE6'   }       0.19817       0.10527     0.32334      250.12
    {'HE7'   }       0.18602       0.10525     0.28267      250.33
    {'JE1'   }       0.21787       0.36942      1.5552         Inf
    {'JE2'   }      0.088423       0.14326     0.53446      659.51
    {'JE3'   }      0.098008       0.18861     0.73237      1613.2
    {'REA1'  }      0.015828      0.066411    0.083633    0.071662
    {'REA2'  }      0.015767      0.065104    0.085021     0.18249
    {'REA3'  }      0.026955      0.078388     0.14976      9.5748
    {'REA4'  }      0.056984      0.072831     0.11414     0.55094
    {'DIS1'  }       0.06049      0.069474    0.094627     0.50523
    {'DIS2'  }      0.014247       0.06466    0.089864    0.062938
    {'DIS3'  }      0.042367       0.06587     0.09173     0.16964
    {'DIS4'  }      0.042051      0.067277      0.0983     0.16183
    {'DIS5'  }      0.033595      0.064165    0.083547    0.068377
    {'TG1'   }      0.028584      0.072761     0.12505      2.6088
    {'AGS'   }      0.034939      0.074652     0.13008      9.3704
    {'WEC1'  }      0.049618      0.074392     0.13312      4.4817
    {'WEC2'  }      0.047218       0.07367     0.14933      4.2072
    {'WEC3'  }      0.042837       0.07503     0.14493      3.8138
    {'HF1'   }        1.0778        99.467      2221.5         Inf
    {'BDT1'  }      0.080257      0.076185     0.23583      2.7487
    {'BDT2'  }        2.0233        10.059      126.44         Inf
    {'MFP'   }      0.029714      0.076232    0.093639    0.075933
    {'UWV'   }      0.088879      0.088331     0.10885      1.0028
    {'IH'    }       0.62393       0.12586     0.24739      288.86
    {'CSE1'  }      0.051537      0.096063     0.22176      102.24
    {'CSE2'  }       0.69479        3.1575      20.655         Inf
    {'EB1'   }       0.01921      0.071743     0.11563     0.90759
    {'EB2'   }      0.017182      0.069954     0.11211     0.90799
    {'EB3'   }      0.017092      0.070808     0.11406      1.2182
    {'EB4'   }      0.057766       0.12968     0.40463      120.76
    {'EB5'   }       0.41112       0.68193      4.1797         Inf
    {'EB6'   }        7.1688        269.53         Inf         Inf
    {'TF1'   }      0.021768      0.072715    0.090289     0.21322
    {'TF2'   }      0.022971      0.068015     0.08977      0.2454
    {'TF3'   }       0.02229      0.071703    0.089702     0.21362
    {'PSM'   }      0.019628      0.071799    0.085939     0.26848
    {'TL'    }       0.11421         10.14         Inf         Inf
    {'CDP'   }        1.3698        61.996      1374.5         Inf
    {'NN1'   }     0.0093114      0.077504     0.15031    0.061514
    {'NN2'   }     0.0086505      0.086895    0.092074    0.061357
    {'NN3'   }      0.011033      0.072068    0.082461    0.067113
    {'NN4'   }      0.015143      0.072042    0.082858    0.069864
    {'NN5'   }      0.023076      0.074356    0.096431     0.33476
    {'NN6'   }      0.027102      0.079382     0.12514     0.96059
    {'NN7'   }      0.027482      0.076597     0.12133     0.96113
    {'NN8'   }      0.016241      0.072468    0.079091    0.063672
    {'NN9'   }      0.026931        0.0711    0.084905    0.089825
    {'NN10'  }      0.038744       0.07687     0.10069     0.55565
    {'NN11'  }      0.082026       0.10463     0.20153        60.4
    {'NN12'  }       0.02017      0.074224    0.088234      0.1278
    {'NN13'  }       0.01815      0.071881    0.092605     0.37041
    {'NN14'  }      0.018332      0.071518    0.090261     0.14175
    {'NN15'  }      0.026882      0.071288    0.087887    0.064997
    {'NN16'  }      0.041634      0.075434    0.088621     0.39144
    {'NN17'  }      0.013669       0.07192    0.080678    0.064377
    {'NN18'  }        86.112           Inf         Inf         Inf
    {'HF2D1' }          9302           Inf         Inf         Inf
    {'HF2D2' }        4339.1           Inf         Inf         Inf
    {'HF2D3' }          8579           Inf         Inf         Inf
    {'HF2D4' }         714.8           Inf         Inf         Inf
    {'HF2D5' }        8670.6           Inf         Inf         Inf
    {'HF2D6' }         689.9           Inf         Inf         Inf
    {'HF2D7' }        9750.3           Inf         Inf         Inf
    {'HF2D8' }        1750.5           Inf         Inf         Inf
    {'HF2D9' }        5972.2           Inf         Inf         Inf
    {'HF2D10'}      0.070926        1.8676      1.4281      1.0676
    {'HF2D11'}      0.015564       0.37285     0.25438     0.29778
    {'HF2D12'}      0.016241       0.29759     0.21215     0.21035
    {'HF2D13'}       0.01764       0.18417     0.14919     0.15091
    {'HF2D14'}      0.027391       0.85806     0.56187     0.51806
    {'HF2D15'}      0.016317       0.12843     0.11362     0.11053
    {'HF2D16'}      0.019536       0.10463    0.098743     0.10089
    {'HF2D17'}      0.015536       0.13288      0.1073     0.10266
    {'HF2D18'}      0.019803       0.12763     0.10389    0.088276
    {'CM1'   }       0.12579        1.0303      1.8974      141.16
    {'CM2'   }       0.90553        4.2086      21.419         Inf
    {'CM3'   }        2.5937        71.627      1907.2         Inf
    {'CM4'   }        19.917        2072.8         Inf         Inf
    {'CM5'   }        92.071           Inf         Inf         Inf
    {'CM6'   }        404.67           Inf         Inf         Inf
    {'TMD'   }      0.023571       0.58392     0.31416     0.34639
    {'FS'    }      0.088428       0.10812     0.75525     0.19589
    {'DLR1'  }      0.043259      0.080158     0.24884      1.6184
    {'DLR2'  }       0.35581       0.49688      3.4813         Inf
    {'DLR3'  }       0.32818       0.48439      3.4235         Inf
    {'ISS1'  }        19.417        2339.7         Inf         Inf
    {'ISS2'  }        18.841        2340.8         Inf         Inf
    {'CBM'   }        38.092           Inf         Inf         Inf
    {'LAH'   }       0.26874        1.0699      7.4225         Inf



The next table provides the optimal values achieved by each solver (larger is better). Note that
the Mosek/SeDuMI solution for EB5 and EB6 is highly inaccurate.

     Problem      customSolver      Mosek         SeDuMi       LMILab  
    __________    ____________    __________    __________    _________

    {'AC1'   }         8.4292         8.4292        8.4292       8.4263
    {'AC2'   }         8.4292         8.4292        8.4292       8.4263
    {'AC3'   }         5.7511          5.751         5.751       5.7495
    {'AC4'   }          407.1          407.1         407.1       406.95
    {'AC5'   }      0.0027711      0.0027696     0.0027696    0.0027692
    {'AC6'   }         116.96         116.96        116.96       116.92
    {'AC7'   }         122.39         122.39        122.39       122.37
    {'AC8'   }         160.77         160.77        160.77       160.72
    {'AC9'   }         198.44         198.44        198.44       198.36
    {'AC10'  }     1.2032e+06      5.802e+05    4.9835e+05          Inf
    {'AC11'  }         108.28         108.28        108.28       108.24
    {'AC12'  }         4.2405          4.239         4.239       4.2373
    {'AC13'  }         4717.8         4717.8        4717.8       4716.3
    {'AC14'  }           6286           6286          6286          Inf
    {'AC15'  }         1.8622          1.862         1.862       1.8618
    {'AC16'  }         1.8622          1.862         1.862       1.8618
    {'AC17'  }         10.951         10.951        10.951       10.946
    {'AC18'  }         856.14         856.14        856.14       855.89
    {'HE1'   }         13.188         13.188        13.188       13.183
    {'HE2'   }         27.672         27.672        27.672       27.662
    {'HE3'   }         8.3866         8.3864        8.3864       8.3833
    {'HE4'   }         38.991         38.991        38.991       38.975
    {'HE5'   }         38.991         38.991        38.991       38.975
    {'HE6'   }          153.6          153.6         153.6       153.55
    {'HE7'   }          153.6          153.6         153.6       153.55
    {'JE1'   }          14752          14752         14752          Inf
    {'JE2'   }         4006.7         4006.6        4006.6       4006.2
    {'JE3'   }         4377.8         4377.7        4377.7       4375.6
    {'REA1'  }         33.549         33.549        33.549       33.537
    {'REA2'  }         33.542         33.542        33.542       33.528
    {'REA3'  }          183.4         183.39        183.39       183.28
    {'REA4'  }       0.056785       0.043371      0.046065          NaN
    {'DIS1'  }         17.246         17.246        17.246       17.242
    {'DIS2'  }         13.514         13.514        13.514        13.51
    {'DIS3'  }         20.306         20.306        20.306       20.302
    {'DIS4'  }          15.43          15.43         15.43       15.426
    {'DIS5'  }      0.0057768      0.0057747     0.0057747    0.0057737
    {'TG1'   }         331.99         331.99        331.99       331.88
    {'AGS'   }         91.715         91.715        91.715       91.669
    {'WEC1'  }         866.58         866.58        866.58       866.19
    {'WEC2'  }         873.75         873.75        873.75       873.24
    {'WEC3'  }         880.24         880.24        880.24       879.97
    {'HF1'   }          67887          67887         67887          Inf
    {'BDT1'  }        0.68685        0.68604       0.68604      0.68579
    {'BDT2'  }         2284.7         2284.7        2284.7          Inf
    {'MFP'   }         2.3356         2.3356        2.3356       2.3346
    {'UWV'   }          75597          75597             1        75582
    {'IH'    }         30.837         30.837        30.837       30.832
    {'CSE1'  }         20.271          20.27         20.27       20.265
    {'CSE2'  }         60.276         29.256        37.922          Inf
    {'EB1'   }          4.332         4.3312        4.3312       4.3309
    {'EB2'   }          4.332         4.3312        4.3312       4.3309
    {'EB3'   }         2.8921         2.8911        2.8911       2.8902
    {'EB4'   }         5.4186         5.3708        5.3708       5.3702
    {'EB5'   }         9.8088     4.8088e-22      2.16e-12          Inf
    {'EB6'   }         38.724     7.4703e-30           Inf          Inf
    {'TF1'   }         5.7969         5.7968        5.7968       5.7939
    {'TF2'   }         5.7969         5.7968        5.7968       5.7939
    {'TF3'   }         5.7969         5.7968        5.7968       5.7939
    {'PSM'   }         70.719         70.719        70.719       70.686
    {'TL'    }              0            NaN           Inf          Inf
    {'CDP'   }          69355          69355         69355          Inf
    {'NN1'   }         7.3099         7.3099        7.3099       7.3085
    {'NN2'   }         1.0308         1.0308        1.0308       1.0306
    {'NN3'   }         1.7708         1.7708        1.7708       1.7702
    {'NN4'   }         10.817         10.817        10.817       10.815
    {'NN5'   }         23.618         23.618        23.618       23.609
    {'NN6'   }         53.164         53.151        53.151       53.119
    {'NN7'   }         53.164         53.151        53.151       53.119
    {'NN8'   }         3.7364         3.7364        3.7364       3.7354
    {'NN9'   }         2.9233         2.9233        2.9233       2.9225
    {'NN10'  }         14.179         14.179        14.179       14.175
    {'NN11'  }         839.49         839.49        839.49       839.19
    {'NN12'  }         6.1863         6.1862        6.1862       6.1831
    {'NN13'  }         44.466         44.465        44.465       44.456
    {'NN14'  }         44.466         44.465        44.465       44.456
    {'NN15'  }         0.2424        0.24195       0.24195       0.2419
    {'NN16'  }         9.4993         9.4992        9.4992       9.4956
    {'NN17'  }         7.0773         7.0773        7.0773        7.076
    {'NN18'  }     1.0011e+06            Inf           Inf          Inf
    {'HF2D1' }       9.63e+06            Inf           Inf          Inf
    {'HF2D2' }      4.107e+07            Inf           Inf          Inf
    {'HF2D3' }     9.5788e+06            Inf           Inf          Inf
    {'HF2D4' }     8.1906e+06            Inf           Inf          Inf
    {'HF2D5' }     9.5752e+06            Inf           Inf          Inf
    {'HF2D6' }      8.183e+06            Inf           Inf          Inf
    {'HF2D7' }     9.5762e+06            Inf           Inf          Inf
    {'HF2D8' }     8.1873e+06            Inf           Inf          Inf
    {'HF2D9' }     9.9922e+05            Inf           Inf          Inf
    {'HF2D10'}         178.57         178.57        178.57       178.54
    {'HF2D11'}         179.23         179.23        179.23       179.15
    {'HF2D12'}         125.76         125.76        125.76       125.73
    {'HF2D13'}         253.87         253.87        253.87       253.77
    {'HF2D14'}         118.31         118.31        118.31       118.27
    {'HF2D15'}         260.24         260.24        260.24       260.11
    {'HF2D16'}         87.677         87.677        87.677       87.643
    {'HF2D17'}          149.9          149.9         149.9       149.85
    {'HF2D18'}         142.36         142.36        142.36       142.32
    {'CM1'   }          4.283          4.283         4.283       4.2812
    {'CM2'   }         10.464         10.464        10.464          Inf
    {'CM3'   }         33.552         33.549        33.549          Inf
    {'CM4'   }         201.37         201.37           Inf          Inf
    {'CM5'   }         1527.6            Inf           Inf          Inf
    {'CM6'   }          12133            Inf           Inf          Inf
    {'TMD'   }         4.1304         4.1303        4.1303       4.1285
    {'FS'    }       0.077008        0.11838      0.091643          NaN
    {'DLR1'  }        0.71914        0.71713       0.71713      0.71666
    {'DLR2'  }         25.379         25.318        25.318          Inf
    {'DLR3'  }         25.379         25.318        25.318          Inf
    {'ISS1'  }         86.834         86.632           Inf          Inf
    {'ISS2'  }         86.834         86.632           Inf          Inf
    {'CBM'   }          12630            Inf           Inf          Inf
    {'LAH'   }         141.33         141.33        141.33          Inf
