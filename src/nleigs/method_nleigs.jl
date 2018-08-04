include("lusolver.jl")
include("spmultiply.jl")
include("discretizepolygon.jl")
include("lejabagby.jl")
include("ratnewtoncoeffs.jl")
include("ratnewtoncoeffsm.jl")

#=
gun_variantP: very accurate: vector V has 11 accurate digits after 100 iterations
gun_variantR1: very accurate: vector V has 11 accurate digits after 100 iterations
gun_variantR2: vector V drifts a bit over time, especially noticeable after ~60 its
gun_variantS: vector V drifts a bit over time, especially noticeable after ~60 its
particle_variantS: vector V drifts more quickly, single digit precision only after ~60 its
particle_variantR2: vector V drifts more quickly, single digit precision only after ~60 its

particle_variantS:
785.781555 seconds (37.02 M allocations: 644.420 GiB, 21.06% gc time)
Found 2 lambdas:
-0.1405984886820147 - 1.817533804824943e-15im
-0.13254142676357075 + 3.230068093242451e-13im

particle_variantR2:
941.252558 seconds (29.10 M allocations: 188.096 GiB, 6.75% gc time)
Found 2 lambdas:
-0.14059848868195352 - 5.4691096223778965e-14im
-0.13254142676311503 - 9.670583417954483e-13im

With reuselu = 2 for variantS, there's still a discrepancy (exactly the same V)
With leja = 0 for variantS, final vector V has 11 accurate digits

Commenting out the two "if leja == 1", V seems to have an even larger discrepancy
So the source of the discrepancy seems to stem from the setup part

julia> workspace(); include("src/nleigs/gun_variantS.jl")
gun_variantS
Starting nleigs
sigma discretization
Rational Newton coefficients
Linearization converged after 35 iterations
 --> freeze linearization
new vector V: size = (12896,), norm = 1.0, sum = -0.16521162667568734 - 0.4412591361800707im
  iteration 1: 0 of 1 < 1.0e-10
new vector V: size = (12896,), norm = 1.0, sum = -0.39708475521655157 + 1.1699737268075565im
  iteration 2: 0 of 0 < 1.0e-10
new vector V: size = (12896,), norm = 1.0, sum = 0.3730074542579329 + 2.2202364594333037im
  iteration 3: 0 of 2 < 1.0e-10
new vector V: size = (12896,), norm = 0.9999999999999999, sum = 0.38142000901911177 - 0.13221205700792815im
  iteration 4: 0 of 1 < 1.0e-10
new vector V: size = (12896,), norm = 0.9999999999999999, sum = 0.5987746981831032 + 0.10704402045719319im
  iteration 5: 0 of 4 < 1.0e-10
new vector V: size = (12896,), norm = 1.0, sum = -1.8292850082266054 + 0.2178943551945447im
  iteration 6: 0 of 5 < 1.0e-10
new vector V: size = (12896,), norm = 1.0, sum = -0.3867354180279049 - 2.3476879366465906im
  iteration 7: 0 of 0 < 1.0e-10
new vector V: size = (12896,), norm = 0.9999999999999999, sum = 1.6774575387042252 - 0.9781846103497509im
  iteration 8: 0 of 2 < 1.0e-10
new vector V: size = (12896,), norm = 0.9999999999999999, sum = -1.9957431343396912 + 0.43268567883549747im
  iteration 9: 0 of 3 < 1.0e-10
new vector V: size = (12896,), norm = 1.0, sum = -0.8710676190417734 - 3.1491850097302887im
  iteration 10: 0 of 4 < 1.0e-10
new vector V: size = (12896,), norm = 1.0, sum = -1.5799152474779754 - 1.5585118887846532im
  iteration 11: 0 of 5 < 1.0e-10
new vector V: size = (12896,), norm = 0.9999999999999999, sum = -0.3814821297398363 + 1.2690750658862602im
  iteration 12: 0 of 4 < 1.0e-10
new vector V: size = (12896,), norm = 0.9999999999999999, sum = -1.0289113160468037 + 3.9573263041926907im
  iteration 13: 0 of 6 < 1.0e-10
new vector V: size = (12896,), norm = 1.0, sum = 1.3714235283754102 + 3.727511609420748im
  iteration 14: 0 of 4 < 1.0e-10
new vector V: size = (12896,), norm = 1.0, sum = -1.160543110389784 - 0.772365326539996im
  iteration 15: 0 of 6 < 1.0e-10
new vector V: size = (12896,), norm = 0.9999999999999999, sum = -1.4385844303839113 - 1.3088980982427945im
  iteration 16: 0 of 7 < 1.0e-10
new vector V: size = (12896,), norm = 1.0, sum = 3.703322212661177 - 3.1719305896847474im
  iteration 17: 0 of 5 < 1.0e-10
new vector V: size = (12896,), norm = 1.0, sum = 3.017327190742389 + 0.8182865636851926im
  iteration 18: 0 of 7 < 1.0e-10
new vector V: size = (12896,), norm = 1.0, sum = 1.5961699048917906 - 0.9408469592608529im
  iteration 19: 0 of 6 < 1.0e-10
new vector V: size = (12896,), norm = 1.0, sum = 1.2233939059824133 - 1.3418538436984544im
  iteration 20: 0 of 7 < 1.0e-10
new vector V: size = (12896,), norm = 0.9999999999999999, sum = -0.840336280313484 + 0.5772831625472992im
  iteration 21: 0 of 7 < 1.0e-10
new vector V: size = (12896,), norm = 1.0, sum = 0.22614191130343264 + 0.7821853212032935im
  iteration 22: 0 of 5 < 1.0e-10
new vector V: size = (12896,), norm = 1.0, sum = -0.03987776353456829 + 1.6557039200701096im
  iteration 23: 0 of 12 < 1.0e-10
new vector V: size = (12896,), norm = 1.0, sum = -0.022844459948222706 + 1.4368096336378218im
  iteration 24: 0 of 11 < 1.0e-10
new vector V: size = (12896,), norm = 1.0000000000000002, sum = -1.512762507106728 - 2.106276989290346im
  iteration 25: 0 of 12 < 1.0e-10
new vector V: size = (12896,), norm = 1.0, sum = 0.3123306865091444 + 0.564125455631111im
  iteration 26: 0 of 14 < 1.0e-10
new vector V: size = (12896,), norm = 1.0, sum = -1.0888716866316202 + 0.45272879663045473im
  iteration 27: 0 of 12 < 1.0e-10
new vector V: size = (12896,), norm = 1.0000000000000002, sum = -0.2691442604077201 + 0.091421540837616im
  iteration 28: 0 of 15 < 1.0e-10
new vector V: size = (12896,), norm = 1.0, sum = -1.5966027018348377 - 1.1179599641584543im
  iteration 29: 0 of 12 < 1.0e-10
new vector V: size = (12896,), norm = 0.9999999999999999, sum = 0.41388251341953375 + 1.1818579315590645im
  iteration 30: 0 of 17 < 1.0e-10
new vector V: size = (12896,), norm = 1.0, sum = 0.1829108445355545 + 0.349279229648053im
  iteration 31: 0 of 15 < 1.0e-10
new vector V: size = (12896,), norm = 0.9999999999999999, sum = -1.3829267814908373 - 1.178454632597906im
  iteration 32: 0 of 14 < 1.0e-10
new vector V: size = (12896,), norm = 1.0, sum = -0.467279843107367 + 0.9770796340227011im
  iteration 33: 0 of 17 < 1.0e-10
new vector V: size = (12896,), norm = 1.0, sum = -1.0905646708473058 - 0.5994135313663044im
  iteration 34: 0 of 15 < 1.0e-10
new vector V: size = (12896,), norm = 1.0, sum = -0.4786916345793831 + 3.685498023360707im
  iteration 35: 0 of 18 < 1.0e-10
new vector V: size = (12896,), norm = 0.9999999999999999, sum = 0.17265141044607812 - 0.4437516303620501im
  iteration 36: 1 of 16 < 1.0e-10
new vector V: size = (12896,), norm = 1.0, sum = -0.25584532419056194 + 2.6709245060237463im
  iteration 37: 1 of 18 < 1.0e-10
new vector V: size = (12896,), norm = 1.0, sum = -1.4905665041275018 - 1.209032353602772im
  iteration 38: 2 of 18 < 1.0e-10
new vector V: size = (12896,), norm = 1.0, sum = 3.1139846928461523 + 1.5474370332835794im
  iteration 39: 2 of 19 < 1.0e-10
new vector V: size = (12896,), norm = 1.0, sum = -1.8365958533424918 + 0.099552506461968im
  iteration 40: 3 of 19 < 1.0e-10
new vector V: size = (12896,), norm = 1.0, sum = -0.9908044044124468 - 0.37594201809220423im
  iteration 41: 1 of 20 < 1.0e-10
new vector V: size = (12896,), norm = 1.0, sum = 0.5776291760328662 + 1.6171308372253974im
  iteration 42: 3 of 19 < 1.0e-10
new vector V: size = (12896,), norm = 1.0, sum = -1.650817087474814 + 1.757491820212293im
  iteration 43: 3 of 20 < 1.0e-10
new vector V: size = (12896,), norm = 1.0, sum = 2.516873463341806 - 1.0691252239545725im
  iteration 44: 4 of 20 < 1.0e-10
new vector V: size = (12896,), norm = 1.0, sum = -3.285328724482927 + 4.202223933941806im
  iteration 45: 4 of 20 < 1.0e-10
new vector V: size = (12896,), norm = 1.0, sum = 0.5358388869221518 - 0.5757060705468127im
  iteration 46: 5 of 21 < 1.0e-10
new vector V: size = (12896,), norm = 1.0, sum = 0.35594318902007305 + 0.23408998936998265im
  iteration 47: 6 of 21 < 1.0e-10
new vector V: size = (12896,), norm = 1.0, sum = -1.0480276337443137 - 1.5838908103465592im
  iteration 48: 9 of 21 < 1.0e-10
new vector V: size = (12896,), norm = 0.9999999999999999, sum = 1.631248237487751 - 1.1930476611931775im
  iteration 49: 11 of 21 < 1.0e-10
new vector V: size = (12896,), norm = 1.0, sum = -3.876825751288304 + 2.0613370501050383im
  iteration 50: 11 of 21 < 1.0e-10
new vector V: size = (12896,), norm = 1.0, sum = 0.8238882683958539 - 0.18042890991615534im
  iteration 51: 11 of 21 < 1.0e-10
new vector V: size = (12896,), norm = 1.0, sum = -0.753857694740055 + 2.7346152880255845im
  iteration 52: 15 of 21 < 1.0e-10
new vector V: size = (12896,), norm = 1.0, sum = 1.255687096027167 + 1.3055223905707336im
  iteration 53: 15 of 21 < 1.0e-10
new vector V: size = (12896,), norm = 1.0, sum = -0.3864992284033473 - 0.9863085874137992im
  iteration 54: 16 of 21 < 1.0e-10
new vector V: size = (12896,), norm = 1.0, sum = 1.2900924341361026 - 1.2520881471017797im
  iteration 55: 16 of 21 < 1.0e-10
new vector V: size = (12896,), norm = 1.0, sum = 1.104759297783835 + 0.7885986365422778im
  iteration 56: 16 of 21 < 1.0e-10
new vector V: size = (12896,), norm = 1.0, sum = -2.2464344160340826 + 0.43289926708715964im
  iteration 57: 16 of 21 < 1.0e-10
new vector V: size = (12896,), norm = 0.9999999999999999, sum = 1.5575882521441053 - 0.9925889017332038im
  iteration 58: 17 of 21 < 1.0e-10
new vector V: size = (12896,), norm = 1.0, sum = -0.7872087162393677 + 1.7273187082443506im
  iteration 59: 17 of 21 < 1.0e-10
new vector V: size = (12896,), norm = 1.0, sum = 1.5564214959655953 - 0.5946919984354854im
  iteration 60: 17 of 21 < 1.0e-10
new vector V: size = (12896,), norm = 1.0, sum = 0.834209843915793 - 0.5854452194521411im
  iteration 61: 17 of 21 < 1.0e-10
new vector V: size = (12896,), norm = 1.0, sum = 0.09427944096889745 - 0.8640834634275937im
  iteration 62: 18 of 21 < 1.0e-10
new vector V: size = (12896,), norm = 1.0, sum = 1.4810647012645293 + 1.4186890555392708im
  iteration 63: 18 of 21 < 1.0e-10
new vector V: size = (12896,), norm = 1.0, sum = -1.0388622072923117 + 1.2429876192760259im
  iteration 64: 19 of 21 < 1.0e-10
new vector V: size = (12896,), norm = 1.0, sum = 0.615107488430894 - 0.6815449582421876im
  iteration 65: 19 of 21 < 1.0e-10
new vector V: size = (12896,), norm = 0.9999999999999999, sum = 0.05072043057994795 - 1.399579284236831im
  iteration 66: 20 of 21 < 1.0e-10
new vector V: size = (12896,), norm = 1.0, sum = 1.286234489479223 + 0.6483612883904784im
  iteration 67: 21 of 21 < 1.0e-10
new vector V: size = (12896,), norm = 0.9999999999999999, sum = 0.1964044102211508 - 0.6972391225419222im
  iteration 68: 21 of 21 < 1.0e-10
new vector V: size = (12896,), norm = 1.0, sum = 0.7458609851760766 - 0.8321983672749984im
  iteration 69: 21 of 21 < 1.0e-10
new vector V: size = (12896,), norm = 1.0, sum = 1.3480684931564784 + 0.16788620331300186im
  iteration 70: 21 of 21 < 1.0e-10

new vector V: size = (12896,), norm = 1.0, sum = 0.615107488430894 - 0.6815449582421876im
new vector V: size = (12896,), norm = 1.0, sum = 0.05072043057994795 - 1.399579284236831im
new vector V: size = (12896,), norm = 1.0, sum = 1.286234489479223 + 0.6483612883904784im
new vector V: size = (12896,), norm = 1.0, sum = 0.1964044102211508 - 0.6972391225419222im
new vector V: size = (12896,), norm = 1.0, sum = 0.7458609851760766 - 0.8321983672749984im
new vector V: size = (12896,), norm = 1.0, sum = 1.3480684931564784 + 0.16788620331300186im

Improved accuracy LU-solver:

new vector V: size = (12896,), norm = 1.0, sum = 0.6253088289271128 - 0.67055241267439im
new vector V: size = (12896,), norm = 1.0, sum = 0.019296820185359928 - 1.3860807195016458im
new vector V: size = (12896,), norm = 1.0, sum = 1.3092181987633762 + 0.6007522684550065im
new vector V: size = (12896,), norm = 1.0, sum = 0.25372344783619344 - 0.6465973735970394im
new vector V: size = (12896,), norm = 1.0, sum = 0.7086570268128205 - 0.751244042471581im
new vector V: size = (12896,), norm = 1.0, sum = 1.14112641420536 + 0.05824430226727742im

High precision V0 + Improved accuracy LU-solver:

new vector V: size = (12896,), norm = 1.0, sum = 0.6130652610061741 - 0.6893661887732254im
new vector V: size = (12896,), norm = 1.0, sum = 0.07068800306178916 - 1.3988670534352587im
new vector V: size = (12896,), norm = 1.0, sum = 1.2682171995504583 + 0.6690848184499923im
new vector V: size = (12896,), norm = 1.0, sum = 0.17643135721400724 - 0.7436463591193164im
new vector V: size = (12896,), norm = 1.0, sum = 0.7857982587112455 - 0.8778716633583774im
new vector V: size = (12896,), norm = 1.0, sum = 1.436559460451974 + 0.3062699452002448im

High precision V0 + Regular LU-solver:

new vector V: size = (12896,), norm = 1.0, sum = 0.6207939109395979 - 0.6882845105377822im
new vector V: size = (12896,), norm = 1.0, sum = 0.06153429305488328 - 1.3854449868729328im
new vector V: size = (12896,), norm = 1.0, sum = 1.2506443067142783 + 0.6391745305233995im
new vector V: size = (12896,), norm = 1.0, sum = 0.2293474041106553 - 0.7444224692266621im
new vector V: size = (12896,), norm = 1.0, sum = 0.8139868560347278 - 0.8088706831693941im
new vector V: size = (12896,), norm = 1.0, sum = 1.2765281238648958 + 0.34241010493088664im

MATLAB:
new vector V: size = (12896,1), norm = 1, sum = 0.621281979143513 + -0.673959376445343i
new vector V: size = (12896,1), norm = 1, sum = 0.0332501670028208 + -1.4044693144825i
new vector V: size = (12896,1), norm = 1, sum = 1.29155225862667 + 0.626507198495382i
new vector V: size = (12896,1), norm = 1, sum = 0.222466890409319 + -0.66169454243909i
new vector V: size = (12896,1), norm = 1, sum = 0.713624580449414 + -0.765114460552119i
new vector V: size = (12896,1), norm = 1, sum = 1.21255281015003 + 0.0754811014799445i

 76.038910 seconds (300.30 M allocations: 20.087 GiB, 17.30% gc time)
Found 21 lambdas:
22345.116783765057 + 0.6449985984624277im
96968.27185274639 + 27532.603459266134im
87004.08355002155 + 28115.999957933218im
44259.41857506462 + 3.5759869516476055im
43857.60089796084 + 20.525532396051673im
54550.13915401979 + 459.51716102748713im
48142.06858696588 + 41.89161304499629im
48788.731987273575 + 6.323940151076352im
109835.02748935619 + 133.7320354120008im
109910.14582387656 + 998.0464371167685im
106301.43146441565 + 86.1611657807596im
106625.99873994477 + 27.03575094242539im
98263.26333958945 + 186.1271754798802im
75402.85310756831 + 4948.348818450018im
77240.79034964913 + 143.90139255653787im
80991.85642218187 + 32.387078392370306im
83158.78304072873 + 458.86690999319086im
88394.7704706994 + 298.72936448365334im
86832.8917008209 + 45.65737695755673im
87407.3563174851 + 35.98153259453339im
87627.51060654216 + 32.130694525689854im
=#

#=
NLEIGS  Find a few eigenvalues and eigenvectors of a NLEP
   lambda = NLEIGS(NLEP,Sigma,[Xi]) returns a vector of eigenvalues of the
   nonlinear eigenvalue problem NLEP inside the target set Sigma. NLEP is a
   structure representing the nonlinear eigenvalue problem as a function
   handle:
     NLEP.Fun:  function handle A = NLEP.Fun(lambda)
     NLEP.n:    size of A
   or as a sum of constant matrices times scalar functions:
     NLEP.B:    array of matrices of the polynomial part
     NLEP.C:    array of matrices of the nonlinear part
     NLEP.L:    array of L-factors of matrices C (optional)
     NLEP.U:    array of U-factors of matrices C (optional)
     NLEP.f:    array of the nonlinear scalar or matrix functions
   Sigma is a vector containing the points of a polygonal target set in the
   complex plane. Xi is a vector containing a discretization of the singularity
   set. In case the singularity set is omitted Xi is set to Inf.

   [X,lambda,res] = NLEIGS(NLEP,Sigma,[Xi]) returns a matrix X of eigenvectors
   and a vector of corresponding eigenvalues of the nonlinear eigenvalue
   problem NLEP inside the target set Sigma. The vector res contains the
   corresponding residuals.

   [X,lambda,res] = NLEIGS(NLEP,Sigma,Xi,options) sets the algorithm's
   parameters to the values in the structure options:
     options.disp:     level of display
                       [ {0} | 1 | 2 ]
     options.maxdgr:   max degree of approximation
                       [ positive integer {100} ]
     options.minit:    min number of iter. after linearization is converged
                       [ positive integer {20} ]
     options.maxit:    max number of total iterations
                       [ positive integer {200} ]
     options.tolres:   tolerance for residual
                       [ positive scalar {1e-10} ]
     options.tollin:   tolerance for convergence of linearization
                       [ positive scalar {1e-11} ]
     options.v0:       starting vector
                       [ vector array {randn(n,1)} ]
     options.funres:   function handle for residual (Lambda[vector], X[matrix])
                       [ @(Lambda,X) {norm(A(lam)*x)} ]
     options.isfunm:   use matrix functions
                       [ boolean {true} ]
     options.static:   static version of nleigs
                       [ boolean {false} ]
     options.leja:     use of Leja-Bagby points
                       [ 0 (no) | {1} (only in expansion phase) | 2 (always) ]
     options.nodes:    prefixed interpolation nodes (only with leja is 0 or 1)
                       [ vector array | {[]} ]
     options.reuselu:  reuse of LU-factorizations of A(sigma)
                       [ 0 (no) | {1} (only after conv. lin) | 2 (always) ]
     options.blksize:  block size for pre-allocation
                       [ positive integer {20} ]

   [X,lambda,res,info] = NLEIGS(NLEP,Sigma,[Xi]) also returns a structure info:
     info.Lam:    matrix of Ritz values in each iteration
     info.Res:    matrix of residuals in each iteraion
     info.sigma:  vector of interpolation nodes
     info.xi:     vector of poles
     info.beta:   vector of scaling parameters
     info.nrmD:   vector of norms of generalized divided differences (in
                  function handle case) or maximum of absolute values of
                  scalar divided differences in each iteration (in matrix
                  function case)

   See also EIG, EIGS, NLEIGSPLOT.

   Reference:
   S. Guettel, R. Van Beeumen, K. Meerbergen, and W. Michiels. NLEIGS: A class
   of fully rational Krylov methods for nonlinear eigenvalue problems. SIAM J.
   Sci. Comput., 36(6), A2842âA2864, 2014.

   Roel Van Beeumen
   April 5, 2016
=#
function nleigs(A, Sigma::Vector{Complex{T}}; Xi = [Inf], options = Dict(), return_info = false) where T

println("Starting nleigs")

# The following variables are used when creating the return values, so put them in scope
# TODO local lam, conv, ...
D = lam = conv = X = res = Lam = Res = sigma = expand = xi = beta = nrmD = maxdgr = kconv = nothing
element_type = eltype(Sigma)

funA,Ahandle,B,BB,pf,C,CC,f,iL,L,LL,U,n,p,q,r,Sigma,leja,nodes,Xi,tollin,
    tolres,maxdgr,minit,maxit,isfunm,static,v0,reuselu,funres,b,computeD,
    resfreq,verbose = checkInputs(A, Sigma, Xi, options)

# Initialization
if computeD
    D = zeros(b+1) # Array{}(b+1) ?
end
if static
    V = zeros(element_type, n, 1)
elseif Ahandle
    V = zeros(element_type, (b+1)*n, b+1)
else
    if r == 0 || b < p
        V = zeros(element_type, (b+1)*n, b+1)
    else
        V = zeros(element_type, p*n+(b-p+1)*r, b+1)
    end
end
H = zeros(element_type, b+1, b)
K = zeros(element_type, b+1, b)
nrmD = Array{Float64}(1)
if return_info
    Lam = zeros(element_type, b, b)
    Res = zeros(b, b)
end

println("sigma discretization")
# Discretization of Sigma --> Gamma & Leja-Bagby points
if leja == 0 # use no leja nodes
    if isempty(nodes)
        error("Interpolation nodes must be provided via 'options[\"nodes\"]' when no Leja-Bagby points (options[\"leja\"] == 0) are used.")
    end
    gamma = discretizepolygon(Sigma)
    max_count = static ? maxit+maxdgr+2 : max(maxit,maxdgr)+2
    sigma = repmat(reshape(nodes, :, 1), ceil(Int, max_count/length(nodes)), 1)
    _,xi,beta = lejabagby(sigma[1:maxdgr+2], Xi, gamma, maxdgr+2, true, p)
elseif leja == 1 # use leja nodes in expansion phase
    if isempty(nodes)
        gamma,nodes = discretizepolygon(Sigma, true)
    else
        gamma = discretizepolygon(Sigma)
        #GammaRealImag = read_sparse_matrix("src/nleigs/gun_gamma.txt")
        #GammaRealImag = full(GammaRealImag)[:]
        #ng = length(GammaRealImag)
        #gamma = GammaRealImag[1:round(Int,ng/2)] + im*GammaRealImag[(round(Int,ng/2)+1):end]
    end
    #@printf("sigma 1569: %s\n", Sigma[1569])
    #@printf("gamma 6098: %s\n", gamma[6098])
    nodes = repmat(reshape(nodes, :, 1), ceil(Int, (maxit+1)/length(nodes)), 1)
    sigma,xi,beta = lejabagby(gamma, Xi, gamma, maxdgr+2, false, p)
    #@printf("gamma: %s\n", norm(gamma))
    #@printf("xi   : %s\n", norm(Xi))
    #@printf("sigma: %s\n", norm(sigma))
#    for i=1:length(sigma)
#        @printf("%d: %s\n", i, sigma[i])
#    end
else # use leja nodes in both phases
    gamma = discretizepolygon(Sigma)
    max_count = static ? maxit+maxdgr+2 : max(maxit,maxdgr)+2
    sigma,xi,beta = lejabagby(gamma, Xi, gamma, max_count, false, p)
end
# verified: gamma, nodes, sigma, xi, beta
xi[maxdgr+2] = NaN # not used
if (Ahandle || !isfunm) && length(sigma) != length(unique(sigma))
    error("All interpolation nodes must be distinct when no matrix " *
        "functions are used for computing the generalized divided differences.")
end

println("Rational Newton coefficients")
# Rational Newton coefficients
range = 1:maxdgr+2
if Ahandle
    D = ratnewtoncoeffs(funA, sigma[range], xi[range], beta[range])
    nrmD[1] = vecnorm(D[1]) # Frobenius norm
else
    # Compute scalar generalized divided differences
    sgddp,sgddq = scgendivdiffs(sigma[range], xi[range], beta[range], pf, f, p, q, maxdgr, isfunm)
    @printf("sgddp: size = %s, norm = %s, sum = %s\n", size(sgddp), norm(sgddp), sum(sum(sgddp)))
    @printf("sgddq: size = %s, norm = %s, sum = %s\n", size(sgddq), norm(sgddq), sum(sum(sgddq)))
    # Construct first generalized divided difference
    if computeD
        D[1] = constructD(0, B, C, L, n, p, q, r, sgddp, sgddq)
    end
    # Norm of first generalized divided difference
    nrmD[1] = maximum(abs.([sgddp[:,1]; sgddq[:,1]]))
end
# verified: sgddp, sgddq, nrmD
if !isfinite(nrmD[1]) # check for NaN
    error("The generalized divided differences must be finite.");
end
lureset()

# Rational Krylov
@printf("v0 A: size = %s, norm = %s, sum = %s\n", size(v0), norm(v0), sum(sum(v0))) # ~15 accurate digits
if reuselu == 2
    v0 = lusolve(funA, sigma[1], v0/norm(v0))
else
    v0 = funA(sigma[1]) \ (v0/norm(v0))
end
#v0 similar, but not exact due to rand
V[1:n,1] = v0/norm(v0)
@printf("v0 B: size = %s, norm = %s, sum = %s\n", size(v0), norm(v0), sum(sum(v0))) # ~14 accurate digits
@printf("V[1:n,1]: size = %s, norm = %s, sum = %s\n", size(V[1:n,1]), norm(V[1:n,1]), sum(sum(V[1:n,1]))) # ~14 accurate digits
expand = true
kconv = Inf
kn = n   # length of vectors in V
l = 0    # number of vectors in V
N = 0    # degree of approximations
kmax = maxit
if static
    kmax += maxdgr
end
k = 1
while k <= kmax
    # allocation
    if l > 0 && (b == 1 || mod(l+1, b) == 1)
        nb = round(Int, 1 + l/b)
        if Ahandle
            V = ensure_size(V, kn+b*n, nb*b+1)
            V[kn+b*n, nb*b+1] = 0
        else
            if expand && computeD
                D[l+b+1] = []
            end
            if expand
                if r == 0 || l + b < p
                    V = ensure_size(V, kn+b*n, nb*b+1)
                    V[kn+b*n, nb*b+1] = 0
                elseif l < p-1
                    V = ensure_size(V, p*n+(nb*b-p+1)*r, nb*b+1)
                    V[p*n+(nb*b-p+1)*r, nb*b+1] = 0
                else # l => p-1
                    V = ensure_size(V, kn+b*r, nb*b+1)
                    V[kn+b*r, nb*b+1] = 0
                end
            else
                V = ensure_size(V, size(V, 1), nb*b+1)
                V[:, nb*b+1] = 0
            end
        end
        # TODO: clean up
        H = ensure_size(H, size(H, 1) + b, size(H, 2) + b)
        H[size(H,1), size(H,2)] = 0
        K = ensure_size(K, size(K, 1) + b, size(K, 2) + b)
        K[size(K,1), size(K,2)] = 0
        if return_info
            Lam = ensure_size(Lam, size(Lam, 1) + b, size(Lam, 2) + b)
            Lam[size(Lam,1), size(Lam,2)] = 0
            Res = ensure_size(Res, size(Res, 1) + b, size(Res, 2) + b)
            Res[size(Res,1), size(Res,2)] = 0
            if expand
                resize!(nrmD, length(nrmD) + b) # TODO: push?
                nrmD[length(nrmD)] = 0
            end
        end
    end

    # set length of vectors
    if expand
        if r == 0 || k < p
            kn += n
        else
            kn += r
        end
    end

    # rational divided differences
    if expand
        if Ahandle
            #
        elseif computeD
            D[k+1] = constructD(k, B, C, L, n, p, q, r, sgddp, sgddq)
        end
        N += 1
    end

    # monitoring norms of divided difference matrices
    if expand
        # TODO: can we pre-allocate nrmD?
        if length(nrmD) < k+1
            resize!(nrmD, k+1)
        end
        if Ahandle
            nrmD[k+1] = vecnorm(D[k+1]) # Frobenius norm
            @printf("nrmD[k+1] (Ahandle) = %s\n", nrmD[k+1])
        else
            # The below can cause out of bounds in sgddp and sgddq if there's
            # no convergence (also happens in MATLAB implementation)
            nrmD[k+1] = maximum(abs.([sgddp[:,k+1]; sgddq[:,k+1]]))
            @printf("nrmD[k+1] = %s\n", nrmD[k+1]) # starts out the same, drifts off slightly over the iterations, likely to be insignificant rounding issues only
        end
        if !isfinite(nrmD[k+1]) # check for NaN
            error("The generalized divided differences must be finite.");
        end
        if n > 1 && k >= 5 && k < kconv
            if sum(nrmD[k-3:k+1]) < 5*tollin
                kconv = k - 1
                if static
                    kmax = maxit + kconv
                end
                expand = false
                if leja == 1
                    # TODO: can we pre-allocate sigma?
                    if length(sigma) < kmax+1
                        resize!(sigma, kmax+1)
                    end
                    sigma[k+1:kmax+1] = nodes[1:kmax-k+1]
                end
                if Ahandle || computeD
                    D = D[1:k]
                end
                xi = xi[1:k]
                beta = beta[1:k]
                nrmD = nrmD[1:k]
                if static
                    if r == 0 || k < p
                        kn -= n
                    else
                        kn -= r
                    end
                    V = ensure_size(V, kn, b+1)
                    V[kn, b+1] = 0
                end
                N -= 1
                if verbose > 0
                    println("Linearization converged after $kconv iterations")
                    println(" --> freeze linearization")
                end
            elseif k == maxdgr+1
                kconv = k
                expand = false
                if leja == 1
                    sigma[k+1:kmax+1] = nodes[1:kmax-k+1]
                end
                if static
                    V = ensure_size(V, kn, b+1)
                    V[kn, b+1] = 0
                end
                N -= 1
                warn("NLEIGS:noconvergence: Linearization not converged after $maxdgr iterations")
                if verbose > 0
                    println(" --> freeze linearization")
                end
            end
        end
    end

    if !static
        l = k
    else
        l = k - N
    end

    # shift-and-invert
    if !static || (static && !expand)
        t = [zeros(l-1,1); 1]  # continuation combination TODO where is this used?
        wc = V[1:kn, l]        # continuation vector
        w = backslash(wc, funA, Ahandle, BB, CC, iL, LL, U, n, p, q, r, reuselu, computeD, sigma, k, D, beta, sgddp, sgddq, N, xi, expand, kconv)
        @printf("kn = %d, l = %d\n", kn, l)
        @printf("wc: size = %s, norm = %s, sum = %s\n", size(wc), norm(wc), sum(sum(wc))) # ~14 accurate digits
        @printf("w: size = %s, norm = %s, sum = %s\n", size(w), norm(w), sum(sum(w))) # ~13 accurate digits
    end

    # orthogonalization
    if !static || (static && !expand)
        normw = norm(w)
        #@printf("normw = %s\n", normw) # differs, first slightly, then more
        h = V[1:kn,1:l]' * w
        w -= V[1:kn,1:l] * h
        H[1:l,l] = h
        eta = 1/sqrt(2)       # reorthogonalization constant
        if norm(w) < eta * normw
            h = V[1:kn,1:l]' * w
            w -= V[1:kn,1:l] * h
            H[1:l,l] += h
        end
        K[1:l,l] = H[1:l,l] * sigma[k+1] + t
    end

    # new vector
    if !static || (static && !expand)
        H[l+1,l] = norm(w)
        #@printf("H[l+1,l] = %s\n", H[l+1,l]) # differ
        K[l+1,l] = H[l+1,l] * sigma[k+1]
        #@printf("K[l+1,l] = %s\n", K[l+1,l]) # differ
        V[1:kn,l+1] = w / H[l+1,l]
        @printf("new vector V: size = %s, norm = %s, sum = %s\n", size(V[1:kn,l+1]), norm(V[1:kn,l+1]), sum(sum(V[1:kn,l+1]))) # ~12 accurate digits
    end

    # Ritz pairs
    if !return_info && ((!expand && k >= N + minit &&
            mod(k-(N+minit), resfreq) == 0) || (k >= kconv + minit &&
            mod(k-(kconv+minit), resfreq) == 0) || k == kmax)
        lambda,S = eig(K[1:l,1:l], H[1:l,1:l])
        # select eigenvalues
        ilam = 1:l
        ilam = ilam[inSigma(lambda, Sigma, tolres)]
        lam = lambda[ilam]
        nblamin = length(lam)
        for i = ilam
            S[:,i] /= norm(H[1:l+1,1:l] * S[:,i])
        end
        X = V[1:n,1:l+1] * (H[1:l+1,1:l] * S[:,ilam])
        for i = 1:size(X,2)
            X[:,i] /= norm(X[:,i])
        end
        # compute residuals
        res = funres(lam, X) # needed for residual: funA, Ahandle, BB, pf, CC, f, p, q
        # check for convergence
        conv = abs.(res) .< tolres
        nbconv = sum(conv)
        if verbose > 0
            iteration = static ? k - N : k
            println("  iteration $iteration: $nbconv of $nblamin < $tolres")
        end
    elseif return_info
        if !static || (static && !expand)
            lambda,S = eig(K[1:l,1:l], H[1:l,1:l])
            # select eigenvalues
            ilam = 1:l
            ilam = ilam[isfinite.(lambda)]
            lam = lambda[ilam]
            lamin = inSigma(lam, Sigma, tolres)
            nblamin = sum(lamin)
            for i = ilam
                S[:,i] /= norm(H[1:l+1,1:l] * S[:,i])
            end
            X = V[1:n,1:l+1] * (H[1:l+1,1:l] * S[:,ilam])
            for i = 1:size(X,2)
                X[:,i] /= norm(X[:,i])
            end
            # compute residuals
            res = zeros(l,1)
            res[isfinite.(lambda)] = NaN
            res[ilam] = funres(lam, X) # needed for residual: funA, Ahandle, BB, pf, CC, f, p, q
            # sort complex numbers by magnitude, then angle
            si = sortperm(lambda, lt = (a,b) -> abs(a) < abs(b) || (abs(a) == abs(b) && angle(a) < angle(b)))
            Res[1:l,l] = res[si]
            Lam[1:l,l] = lambda[si]
            # check for convergence
            res = res[ilam]
            conv = (abs.(res) .< tolres) .& lamin
            nbconv = sum(conv)
            if verbose > 0
                iteration = static ? k - N : k
                println("  iteration $iteration: $nbconv of $nblamin < $tolres")
            end
        end
    end

    # stopping
    if ((!expand && k >= N + minit) || k >= kconv + minit) && nblamin == nbconv
        break
    end

    # increment k
    k += 1
end
lureset()

info = Dict()

if return_info
    Lam = Lam[1:l,1:l]
    Res = Res[1:l,1:l]
    sigma = sigma[1:k]
    if expand
        xi = xi[1:k]
        beta = beta[1:k]
        nrmD = nrmD[1:k]
        warn("NLEIGS:noconvergence: Linearization not converged after $maxdgr iterations")
    end
    info = Dict(
        "Lam" => Lam,
        "Res" => Res,
        "sigma" => sigma,
        "xi" => xi,
        "beta" => beta,
        "nrmD" => nrmD,
        "kconv" => kconv)
end

return X[:,conv], lam[conv], res[conv], info
end
# ------------------------------------------------------------------------------
# Nested functions
# ------------------------------------------------------------------------------

# checkInputs: error checks the inputs to NLEP and also derives some variables
# '''''''''''' from them:
#
#   funA       function handle for A(lambda)
#   Ahandle    is true if NLEP is given as function handle
#   B          cell array {B0,B1,...,Bp}
#   BB         matrix [B0;B1;...;Bp]
#   pf         cell array {x^0,x^1,...,x^p}
#   C          cell array {C1,C2,...,Cq}
#   CC         matrix [C1;C2;...;Cq]
#   f          cell array {f1,f2,...,fq}
#   iL         vector with indices of L-factors
#   L          cell array {L1,L2,...,Lq}
#   LL         matrix [L1,L2,...,Lq]
#   U          matrix [U1,U2,...,Uq]
#   n          dimension of NLEP
#   p          order of polynomial part (length of B = p + 1)
#   q          length of f and C
#   r          sum of ranks of low rank matrices in C
#   Sigma      vector with target set (polygon)
#   leja       0: no leja, 1: leja in expansion phase, 2: always leja
#   nodes      [] or given nodes
#   Xi         vector with singularity set (discretized)
#   computeD   is true if generalized divided differences are explicitly used
#   tollin     tolerance for convergence of linearization
#   tolres     tolerance for residual
#   maxdgr     maximum degree of approximation
#   minit      minimum number of iterations after linearization is converged
#   maxit      maximum number of iterations
#   isfunm     is true if f are matrix functions
#   static     is true if static version is used
#   v0         starting vector
#   reuselu    positive integer for reuse of LU-factorizations of A(sigma)
#   funres     function handle for residual R(Lambda,X)
#   b          block size for pre-allocation
#   verbose    level of display [ {0} | 1 | 2 ]
    function checkInputs(A, Sigma, Xi, options)
        ## initialize
        B = []
        BB = []
        C = []
        CC = []
        pf = []
        f = []
        iL = 0
        L = []
        LL = []
        U = []
        p = -1
        q = 0
        r = 0

        #temp
        funA = []
        residual = []
        computeD = []
        resfreq = []
        #end temp

        ## process the input A
        if !isa(A, Dict)
            error("The input 'A' must be a Dict representing the NLEP.")
        end
        if haskey(A, "Fun")
            Ahandle = true
            # function handle for A(lambda)
            funA = A["Fun"]
        else
            Ahandle = false
            # polynomial part B
            B = get(A, "B", "missing")
            if isempty(B)
                p = -1
            elseif isa(B, Number)
                p = 0;
                B = [B]
            elseif isa(B, Vector)
                p = length(B) - 1
            else
                error("The constant matrices of the polynomial part 'B' must be a 1 dimensional array.")
            end
            # nonlinear part C
            if haskey(A, "C")
                C = A["C"]
                if isempty(C)
                    qC = 0
                elseif isa(C, Number)
                    qC = 1
                    C = [C]
                elseif isa(C, Vector)
                    qC = length(C)
                else
                    error("The constant matrices of the nonlinear part 'C' must be a 1 dimensional array.")
                end
                f = A["f"]
            end
            # L and U factors of the low rank nonlinear part C
            if haskey(A, "L")
                L = A["L"]
                if isa(L, Number)
                    qL = 1
                    L = [L]
                elseif isa(L, Vector)
                    qL = length(L)
                else
                    error("The 'L-factor matrices of the nonlinear part 'C' must be a 1 dimensional array.")
                end
                U = A["U"]
                if isa(U, Number)
                    qU = 1
                    U = [U]
                elseif isa(U, Vector)
                    qU = length(U)
                else
                    error("The 'U'-factor matrices of the nonlinear part 'C' must be a 1 dimensional array.")
                end
                if qL != qU
                    error("The number of 'L'- and 'U'-factors must be equal.")
                end
                if !haskey(A, "C")
                    qC = qL
                    C = zeros(qC)
                    for ii = 1:qC
                        C[ii] = L[ii]*U[ii]'
                    end
                elseif qC != qL
                    error("The number of 'L'- and 'U'-factors must be equal to the number of constant matrices of the nonlinear part 'C'.")
                end
                U = hcat(U...) # input = 81 x 16281 x 2; output = 16281 x 162
                r = size(U, 2)
                iL = zeros(Int, r, 1)
                c = 0
                for ii = 1:qL
                    ri = size(L[ii], 2)
                    iL[c+1:c+ri] = ii
                    c += ri
                end
                LL = hcat(L...)
            end
            # nonlinear function f
            if haskey(A, "C") || haskey(A, "L")
                f = A["f"]
                if isempty(f)
                    qf = 0
                elseif !isempty(methods(f))
                    qf = 1
                    f = [f]
                elseif isa(f, Vector)
                    qf = length(f)
                else
                    error("The scalar functions 'f' of the nonlinear part 'C' must be a 1 dimensional array.")
                end
                if qC == qf
                    q = qC
                else
                    error("The number of constant matrices 'C' and the number of scalar functions 'f' must be equal.")
                end
            end
        end

        ## process the input Xi
        if isempty(Xi)
            Xi = [Inf]
        elseif !isa(Xi, Vector) || !(eltype(Xi) <: Number)
            error("The poles set 'Xi' must be a numeric vector.")
        end

        ## set n and funA
        if Ahandle
            # dimension of A(lambda)
            n = get(A, "n", "missing")
            if !isscalar(n) || !isreal(n) || n < 0 || !isfinite(n)
                error("Size of problem 'n' must be a positive integer.")
            end
            if issparse(n)
                n = full(n)
            end
            if round(n) != n
                warn("WarnTests:convertTest: Size of problem 'n' must be a positive integer. Rounding input size.")
                n = round(n)
            end
        else
            if p >= 0
                n = size(B[1], 1)
                pf = monomials(p)
            elseif q > 0
                n = size(C[1], 1)
            else
                error("'B' or 'C' must be non-empty.")
            end
            if q == 0
                if issparse(B[1])
                    Amat = hcat(map(x -> x[:], B)...) # in: n * x * y; out: (x*y) * n
                else
                    Amat = cat(3, B...) # in: n * x * y; out: x * y * n
                end
                F = pf
            elseif p < 0
                if issparse(C[1])
                    Amat = hcat(map(x -> x[:], C)...) # in: n * x * y, out: (x*y) * n
                else
                    Amat = cat(3, C...) # in: n * x * y; out: x * y * n
                end
                F = f[:]
            else
                if issparse(B[1])
                    Amat = hcat(map(x -> x[:], B)..., map(x -> x[:], C)...)
                else
                    Amat = cat(3, B..., C...)
                end
                F = [pf[:]; f[:]]
            end
            if issparse(Amat)
                funA = lambda -> sparse(reshape(spmultiply(Amat, map(f -> f(lambda), F)), n, n))
            else
                funA = lambda -> sum(Amat[:,:,i] * F[i](lambda) for i=1:length(F))
            end
            BB = vcat(B...)
            CC = vcat(C...)
        end

        ## process the input options
        if !isa(options, Dict)
            error("The input argument 'options' must be a Dict.");
        end

        # extract options, with default values if missing
        verbose = get(options, "disp", 0)
        maxdgr = get(options, "maxdgr", 100)
        minit = get(options, "minit", 20)
        maxit = get(options, "maxit", 200)
        tolres = get(options, "tolres", 1e-10)
        tollin = get(options, "tollin", max(tolres/10, 100*eps()))
        v0 = get(options, "v0", randn(n))
        funres = get(options, "funres", residual)
        isfunm = get(options, "isfunm", true)
        static = get(options, "static", false)
        leja = get(options, "leja", 1)
        nodes = get(options, "nodes", [])
        reuselu = get(options, "reuselu", 1)
        b = get(options, "blksize", 20)

        # extra defaults
        computeD = (n <= 400)
        resfreq = 5

        if n == 1
            maxdgr = maxit + 1;
        end

        return funA,Ahandle,B,BB,pf,C,CC,f,iL,L,LL,U,n,p,q,r,Sigma,leja,nodes,
                Xi,tollin,tolres,maxdgr,minit,maxit,isfunm,static,v0,reuselu,
                funres,b,computeD,resfreq,verbose
    end # checkInputs

# ------------------------------------------------------------------------------

# scgendivdiffs: compute scalar generalized divided differences
#   sigma   discretization of target set
#   xi      discretization of singularity set
#   beta    scaling factors
function scgendivdiffs(sigma, xi, beta, pf, f, p, q, maxdgr, isfunm)
    ## Compute scalar generalized divided differences of polynomials part
    sgddp = complex(zeros(p+1, maxdgr+2))
    for ii = 1:p+1
        if isfunm
            sgddp[ii,:] = ratnewtoncoeffsm(pf[ii], sigma, xi, beta)
        else
            sgddp[ii,:] = ratnewtoncoeffs(pf[ii], sigma, xi, beta)
        end
    end
    ## Compute scalar generalized divided differences of nonlinear part
    sgddq = complex(zeros(q, maxdgr+2))
    for ii = 1:q
        if isfunm
            sgddq[ii,:] = ratnewtoncoeffsm(f[ii], sigma, xi, beta)
            #@printf("ii=%d, sgddq=%s\n", ii, norm(sgddq[ii,:]))
            #@printf("sigma=%s\n", norm(sigma))
            #@printf("beta=%s\n", norm(beta))
        else
            sgddq[ii,:] = ratnewtoncoeffs(f[ii], sigma, xi, beta)
        end
    end

    return sgddp,sgddq
end # scgendivdiffs

# ------------------------------------------------------------------------------

# constructD: Construct generalized divided difference
#   nb  number
    function constructD(nb, B, C, L, n, p, q, r, sgddp, sgddq)
        if r == 0 || nb <= p
            D = spzeros(n, n)
            for ii = 1:p+1
                D += sgddp[ii,nb+1] * B[ii]
            end
            for ii = 1:q
                D += sgddq[ii,nb+1] * C[ii]
            end
        else
            D = []
            for ii = 1:q
                d = sgddq[ii,nb+1] * L[ii]
                D = ii == 1 ? d : hcat(D, d)
            end
        end
        return D
    end # constructD

# ------------------------------------------------------------------------------

# backslash: Backslash or left matrix divide
#   wc       continuation vector
    function backslash(wc, funA, Ahandle, BB, CC, iL, LL, U, n, p, q, r, reuselu, computeD, sigma, k, D, beta, sgddp, sgddq, N, xi, expand, kconv)
        shift = sigma[k+1]
        @printf("shift = %s; k+1 = %d\n", shift, k+1) # same

        ## construction of B*wc
        Bw = zeros(eltype(wc), size(wc))
        # first block (if low rank)
        if r > 0
            i0b = (p-1)*n + 1
            i0e = p*n
            @printf("i0b = %d, i0e = %d\n", i0b, i0e) # same
            if Ahandle || computeD
                Bw[1:n] = -D[p+1] * wc[i0b:i0e] / beta[p+1]
            else
                Bw[1:n] =
                    -sum(reshape(BB * wc[i0b:i0e], n, p+1) .* sgddp[:, p+1].', 2) / beta[p+1] -
                     sum(reshape(CC * wc[i0b:i0e], n, q) .* sgddq[:, p+1].', 2) / beta[p+1]
                @printf("Bw[1:n]: size = %s, norm = %s, sum = %s\n", size(Bw[1:n]), norm(Bw[1:n]), sum(sum(Bw[1:n]))) # same
            end
        end
        # other blocks
        i0b = 1
        i0e = n
        for ii = 1:N
            # range of block i+1
            i1b = i0e + 1
            if r == 0 || ii < p
                i1e = i0e + n
            else
                i1e = i0e + r
            end
            # compute block i+1
            if r == 0 || ii != p
                Bw[i1b:i1e] = wc[i0b:i0e] + beta[ii+1]/xi[ii]*wc[i1b:i1e]
            else
                Bw[i1b:i1e] = U' * wc[i0b:i0e] + beta[ii+1]/xi[ii]*wc[i1b:i1e]
#                @printf("ii = %d; U': sum = %s\n", ii, sum(U'))
#                @printf("ii = %d; wc[i0b:i0e]: norm = %s, sum = %s\n", ii, norm(wc[i0b:i0e]), sum(wc[i0b:i0e]))
#                @printf("ii = %d; U' * wc[i0b:i0e]: norm = %s, sum = %s\n", ii, norm(U' * wc[i0b:i0e]), sum(U' * wc[i0b:i0e]))
#                UUj = real.(U' * wc[i0b:i0e])
#                for ui in eachindex(UUj)
#                    @printf("Uwc %d = %s, abs diff = %s, rel diff = %s\n", ui, UUj[ui], UUj[ui] - UU[ui], (UUj[ui] - UU[ui])/UUj[ui])
#                end
#                @printf("wc[i1b:i1e]: size = %s, norm = %s, sum = %s\n", size(wc[i1b:i1e]), norm(wc[i1b:i1e]), sum(sum(wc[i1b:i1e])))
#                @printf("ii = %d; i1b=%d, i1e=%d, Bw[i1b:i1e]: size = %s, norm = %s, sum = %s\n", ii, i1b, i1e, size(Bw[i1b:i1e]), norm(Bw[i1b:i1e]), sum(sum(Bw[i1b:i1e])))
            end
            # range of next block i
            i0b = i1b
            i0e = i1e
        end

        ## construction of z0
        z = copy(Bw)
        i1b = n + 1
        if r == 0 || p > 1
            i1e = 2*n
        else
            i1e = n + r
        end
        nu = beta[2]*(1 - shift/xi[1])
        #@printf("A z[i1b:i1e]: size = %s, norm = %s, sum = %s\n", size(z[i1b:i1e]), norm(z[i1b:i1e]), sum(sum(z[i1b:i1e])))
        z[i1b:i1e] = 1/nu * z[i1b:i1e]
        #@printf("B z[i1b:i1e]: size = %s, norm = %s, sum = %s\n", size(z[i1b:i1e]), norm(z[i1b:i1e]), sum(sum(z[i1b:i1e])))
        for ii = 1:N
            # range of block i+2
            i2b = i1e + 1
            if r == 0 || ii < p-1
                i2e = i1e + n
            else
                i2e = i1e + r
            end
            # add extra term to z0
            if Ahandle || computeD
                if r == 0 || ii != p
                    z[1:n] -= D[ii+1] * z[i1b:i1e]
                end
            else
                if r == 0 || ii < p
                    if p >= 0
                        z[1:n] -= sum(reshape(BB * z[i1b:i1e], n, p+1) .* sgddp[:,ii+1].', 2)
                    end
                    if q > 0
                        z[1:n] -= sum(reshape(CC * z[i1b:i1e], n, q) .* sgddq[:,ii+1].', 2)
                    end
                elseif ii > p
                    dd = sgddq[:,ii+1]
                    #@printf("A z[1:n]: size = %s, norm = %s, sum = %s\n", size(z[1:n]), norm(z[1:n]), sum(sum(z[1:n])))
                    #@printf("LL: size = %s, sum = %s\n", size(LL), sum(sum(LL)))
                    # this differs on it2:
                    #@printf("i1b=%d, i1e=%d, z[i1b:i1e]: size = %s, norm = %s, sum = %s\n", i1b, i1e, size(z[i1b:i1e]), norm(z[i1b:i1e]), sum(sum(z[i1b:i1e])))
                    #@printf("dd[iL]: size = %s, norm = %s, sum = %s\n", size(dd[iL]), norm(dd[iL]), sum(sum(dd[iL])))
                    z[1:n] -= LL*(z[i1b:i1e] .* dd[iL])
                    #@printf("B z[1:n]: size = %s, norm = %s, sum = %s\n", size(z[1:n]), norm(z[1:n]), sum(sum(z[1:n])))
                end
            end
            # update block i+2
            if ii < N
                mu = shift - sigma[ii+1]
                nu = beta[ii+2] * (1 - shift/xi[ii+1])
                if r == 0 || ii != p-1
                    #@printf("A mu=%s, nu=%s\n", mu, nu)
                    #@printf("A i2b=%d, i2e=%d, z[i2b:i2e]: size = %s, norm = %s, sum = %s\n", i2b, i2e, size(z[i2b:i2e]), norm(z[i2b:i2e]), sum(sum(z[i2b:i2e])))
                    #@printf("A i1b=%d, i1e=%d, z[i1b:i1e]: size = %s, norm = %s, sum = %s\n", i1b, i1e, size(z[i1b:i1e]), norm(z[i1b:i1e]), sum(sum(z[i1b:i1e])))
                    z[i2b:i2e] = 1/nu * z[i2b:i2e] + mu/nu * z[i1b:i1e]
                    #@printf("B i2b=%d, i2e=%d, z[i2b:i2e]: size = %s, norm = %s, sum = %s\n", i2b, i2e, size(z[i2b:i2e]), norm(z[i2b:i2e]), sum(sum(z[i2b:i2e])))
                else # i == p-1
                    z[i2b:i2e] = 1/nu * z[i2b:i2e] + mu/nu * U'*z[i1b:i1e]
                end
            end
            # range of next block i+1
            i1b = i2b
            i1e = i2e
        end

        ## solving Alam x0 = z0
        w = zeros(eltype(wc), size(wc))
        if ((!expand || k > kconv) && reuselu == 1) || reuselu == 2
            w[1:n] = lusolve(funA, shift, z[1:n]/beta[1])
            # the below differs on the 2nd iteration (1st is same) (note: lusolve above is a cache miss)
            #@printf("z[1:n]: size = %s, norm = %s, sum = %s\n", size(z[1:n]), norm(z[1:n]), sum(sum(z[1:n])))
            #@printf("beta[1]: size = %s, norm = %s, sum = %s\n", size(beta[1]), norm(beta[1]), sum(sum(beta[1])))
            #@printf("z[1:n]/beta[1]: size = %s, norm = %s, sum = %s\n", size(z[1:n]/beta[1]), norm(z[1:n]/beta[1]), sum(sum(z[1:n]/beta[1])))
            #@printf("w[1:n]: size = %s, norm = %s, sum = %s\n", size(w[1:n]), norm(w[1:n]), sum(sum(w[1:n])))
        else
            w[1:n] = funA(shift) \ (z[1:n]/beta[1])
        end

        ## substitutions x[i+1] = mu/nu*x[i] + 1/nu*Bw[i+1]
        i0b = 1
        i0e = n
        for ii = 1:N
            # range of block i+1
            i1b = i0e + 1
            if r == 0 || ii < p
                i1e = i0e + n
            else
                i1e = i0e + r
            end
            # compute block i+1
            mu = shift - sigma[ii]
            nu = beta[ii+1] * (1 - shift/xi[ii])
            if r == 0 || ii != p
                w[i1b:i1e] = mu/nu * w[i0b:i0e] + 1/nu * Bw[i1b:i1e]
            else
                w[i1b:i1e] = mu/nu * U'*w[i0b:i0e] + 1/nu * Bw[i1b:i1e]
                #@printf("i0b = %d, i0e = %d, i1b = %d, i1e = %d, w[i1b:i1e]: norm = %s\n", i0b, i0e, i1b, i1e, norm(w[i1b:i1e]))
                #@printf("mu = %s, nu = %s\n", mu, nu)
                #@printf("U'*w[i0b:i0e] norm = %s\n", norm(U'*w[i0b:i0e]))
                #@printf("Bw norm = %s\n", norm(Bw[i1b:i1e])) # much smaller
            end
            # range of next block i
            i0b = i1b
            i0e = i1e
        end

        return w
    end # backslash

# ------------------------------------------------------------------------------

# inSigma: True for points inside Sigma
#   z      (complex) points
    function inSigma(z, Sigma, tolres)
        if length(Sigma) == 2 && isreal(Sigma)
            realSigma = real([Sigma[1]; Sigma[1]; Sigma[2]; Sigma[2]]) # note: sigma may be of complex type with 0 complex part
            imagSigma = [-tolres; tolres; tolres; -tolres]
        else
            realSigma = real(Sigma)
            imagSigma = imag(Sigma)
        end
        return map(p -> inpolygon(real(p), imag(p), realSigma, imagSigma), z)
    end
#=
# ------------------------------------------------------------------------------

# residual: Residual of the NLEP for given eigenvalues and eigenvectors
#   Lambda  eigenvalues
#   X       eigenvectors (normalized)
    function residual(Lambda, X, funA, Ahandle, BB, pf, CC, f, p, q)
        R = zeros(length(Lambda), 1)
        if Ahandle
            for ii = 1:length(Lambda)
                R[ii] = norm(funA(Lambda[ii]) * X[:,ii])
            end
        else
            if p < 0
                BCX = reshape(CC*X, [], q, size(X,2))
            else
                BCX = reshape([BB*X; CC*X], [], p+q+1, size(X,2))
            end
            FL = cell2mat(cellfun(@(x) arrayfun(x,Lambda),[pf[:]; f[:]]', 'UniformOutput',0))
            RR = sum(bsxfun(@times, BCX, reshape(FL.', 1, size(BCX, 2), [])), 2)
            for ii = 1:length(Lambda)
                R[ii] = norm(RR[:, :, ii])
            end
        end
        return R
    end # residual
end # nleigs
=#
function monomials(p)
    return map(i -> x -> x^(i-1), 1:p+1)
end

function ensure_size(A, rows, cols)
    r = max(size(A, 1), rows)
    c = max(size(A, 2), cols)
    if r > size(A, 1) || c > size(A, 2)
        expanded = zeros(eltype(A), r, c)
        expanded[1:size(A, 1), 1:size(A, 2)] = A
        A = expanded
    end
    return A
end
