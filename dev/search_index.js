var documenterSearchIndex = {"docs":
[{"location":"public_api/heuristics/#Heuristic-algorithms","page":"Heuristic solvers","title":"Heuristic algorithms","text":"","category":"section"},{"location":"public_api/heuristics/","page":"Heuristic solvers","title":"Heuristic solvers","text":"optimalportfolio_greedy\noptimalportfolio_simulatedannealing","category":"page"},{"location":"public_api/heuristics/#OptimalApplication.optimalportfolio_greedy","page":"Heuristic solvers","title":"OptimalApplication.optimalportfolio_greedy","text":"optimalportfolio_greedy(mkt::VariedCostsMarket) -> X, v\n\nUse a greedy heuristic that adds schools in decreasing order of mkt.ft ./ mkt.g to compute a heuristically optimal portfolio for the VariedCostsMarket defined by mkt.\n\n\n\n\n\n","category":"function"},{"location":"public_api/heuristics/#OptimalApplication.optimalportfolio_simulatedannealing","page":"Heuristic solvers","title":"OptimalApplication.optimalportfolio_simulatedannealing","text":"optimalportfolio_simulatedannealing(\n    mkt::VariedCostsMarket;\n    temp::Float64=0.25,\n    red::Float64=0.0625,\n    nit::Integer=500,\n    verbose::Bool=false\n) -> X, v\n\nUse a simulated annealing procedure to compute a heuristically optimal portfolio for the VariedCostsMarket defined by mkt. temp is the initial temperature, red is the reduction factor, and nit is the total number of iterations.\n\n\n\n\n\n","category":"function"},{"location":"examples/reach_safety/#Reach-schools-and-safety-schools","page":"Reach & safety schools","title":"Reach schools and safety schools","text":"","category":"section"},{"location":"examples/reach_safety/","page":"Reach & safety schools","title":"Reach & safety schools","text":"Let's compute the optimal portfolio for a VariedCostsMarket, plot the results, and see if we can see any patterns.","category":"page"},{"location":"examples/reach_safety/","page":"Reach & safety schools","title":"Reach & safety schools","text":"First we will read in our packages and generate a random market.","category":"page"},{"location":"examples/reach_safety/","page":"Reach & safety schools","title":"Reach & safety schools","text":"using UnicodePlots, Crayons\nusing OptimalApplication\nmkt = VariedCostsMarket(30)\nscatterplot(\n    mkt.f,\n    mkt.t,\n    marker=[k == 10 ? \"T\" : \"$k\" for k in mkt.g],\n    color=crayon\"yellow\",\n    xlabel=\"f\",\n    ylabel=\"t\",\n    width=60,\n    height=25\n)","category":"page"},{"location":"examples/reach_safety/","page":"Reach & safety schools","title":"Reach & safety schools","text":"As the plot shows, calling VariedCostsMarket(m) generates a random problem in which f and t are negatively correlated. This makes the problem more interesting. By default, the budget mkt.H is chosen so that about half the schools end up in the optimal portfolio. (The numbers in the plot indicate the application cost, where T means \"ten.\")","category":"page"},{"location":"examples/reach_safety/","page":"Reach & safety schools","title":"Reach & safety schools","text":"We can find the optimal solution and objective value using optimalportfolio_dynamicprogram.","category":"page"},{"location":"examples/reach_safety/","page":"Reach & safety schools","title":"Reach & safety schools","text":"X, v = optimalportfolio_dynamicprogram(mkt)","category":"page"},{"location":"examples/reach_safety/","page":"Reach & safety schools","title":"Reach & safety schools","text":"Now let's plot the results:","category":"page"},{"location":"examples/reach_safety/","page":"Reach & safety schools","title":"Reach & safety schools","text":"X_comp = setdiff(1:mkt.m, X)\npl = scatterplot(\n    mkt.f[X_comp],\n    mkt.t[X_comp],\n    marker=[k == 10 ? \"T\" : \"$k\" for k in mkt.g[X_comp]],\n    color=crayon\"red\",\n    name=\"don't apply\",\n    xlabel=\"f\",\n    ylabel=\"t\",\n    width=60,\n    height=25 \n)\nscatterplot!(\n    pl,\n    mkt.f[X],\n    mkt.t[X],\n    marker=[k == 10 ? \"T\" : \"$k\" for k in mkt.g[X]],\n    color=crayon\"green\",\n    name=\"apply\"\n)","category":"page"},{"location":"examples/reach_safety/","page":"Reach & safety schools","title":"Reach & safety schools","text":"The student should apply to the schools shown in green. For a large budget like this, schools with high utility and low admissions probability (\"reach schools\") tend to be preferable to schools with low utility and high admissions probability (\"safety schools\") because the number of schools in the portfolio (and therefore the odds of getting at least one acceptance) is high.","category":"page"},{"location":"examples/reach_safety/","page":"Reach & safety schools","title":"Reach & safety schools","text":"note: Note\nIn creating the plot above, mkt.f[X] and so on retrieve the data in the correct permutation because the random constructor ensures that mkt.perm = 1:mkt.m. For a general input data, you need to use mkt.f[invperm(mkt.perm)[X]], or simply save f and t beforehand. ","category":"page"},{"location":"public_api/approximation/#Approximation-algorithms","page":"Approximation methods","title":"Approximation algorithms","text":"","category":"section"},{"location":"public_api/approximation/","page":"Approximation methods","title":"Approximation methods","text":"optimalportfolio_fptas","category":"page"},{"location":"public_api/approximation/#OptimalApplication.optimalportfolio_fptas","page":"Approximation methods","title":"OptimalApplication.optimalportfolio_fptas","text":"optimalportfolio_fptas(mkt, ε; verbose=false) -> X, v\n\nUse the fully polynomial-time approximation scheme to produce a 1-ε-optimal portfolio X, with valuation v, for the VariedCostsMarket defined by mkt.\n\n\n\n\n\n","category":"function"},{"location":"public_api/fast_solvers/#Fast-solvers","page":"Fast, exact solvers","title":"Fast solvers","text":"","category":"section"},{"location":"public_api/fast_solvers/","page":"Fast, exact solvers","title":"Fast, exact solvers","text":"applicationorder_list\napplicationorder_heap\noptimalportfolio_dynamicprogram","category":"page"},{"location":"public_api/fast_solvers/#OptimalApplication.applicationorder_list","page":"Fast, exact solvers","title":"OptimalApplication.applicationorder_list","text":"applicationorder_list(mkt::SameCostsMarket) -> X, V\n\nProduce the optimal application order X and associated valuations V  for the SameCostsMarket defined by mkt. Uses a list data structure; typically faster than the equivalent applicationorder_heap.\n\nAll SameCostsMarkets satisfy a nestedness property, meaning that the optimal portfolio when mkt.h = h is given by the first h entries of the optimal portfolio when mkt.h = mkt.m. For example:\n\njulia> mkt = SameCostsMarket([0.2, 0.5, 0.1, 0.6, 0.1], [1, 4, 9, 1, 8], 5);\n\njulia> x, v = applicationorder_list(mkt)\n([2, 3, 5, 4, 1], [2.0, 2.7, 3.24, 3.483, 3.5154])\n\njulia> x[1:4], v[4] # optimal portfolio and valuation for h = 4 \n([2, 3, 5, 4], 3.483)\n\n\n\n\n\n","category":"function"},{"location":"public_api/fast_solvers/#OptimalApplication.applicationorder_heap","page":"Fast, exact solvers","title":"OptimalApplication.applicationorder_heap","text":"applicationorder_list(mkt::SameCostsMarket) -> X, V\n\nProduce the optimal application order X and associated valuations V  for the SameCostsMarket defined by mkt. Uses a heap data structure; typically the equivalent applicationorder_list is faster.\n\nAll SameCostsMarkets satisfy a nestedness property, meaning that the optimal portfolio when mkt.h = h is given by the first h entries of the optimal portfolio when mkt.h = mkt.m. For example:\n\njulia> mkt = SameCostsMarket([0.2, 0.5, 0.1, 0.6, 0.1], [1, 4, 9, 1, 8], 5);\n\njulia> x, v = applicationorder_heap(mkt)\n([2, 3, 5, 4, 1], [2.0, 2.7, 3.24, 3.483, 3.5154])\n\njulia> x[1:4], v[4] # optimal portfolio and valuation for h = 4 \n([2, 3, 5, 4], 3.483)\n\n\n\n\n\n","category":"function"},{"location":"public_api/fast_solvers/#OptimalApplication.optimalportfolio_dynamicprogram","page":"Fast, exact solvers","title":"OptimalApplication.optimalportfolio_dynamicprogram","text":"optimalportfolio_dynamicprogram(mkt::VariedCostsMarket) -> X, v\n\nUse the dynamic program on application costs to produce the optimal portfolio X and associated value v for the VariedCostsMarket defined by mkt. \n\njulia> mkt = VariedCostsMarket([0.2, 0.5, 0.1, 0.6, 0.1], [1, 4, 9, 1, 8], [2, 4, 2, 5, 1], 8);\n\njulia> optimalportfolio_dynamicprogram(mkt)\n([3, 5, 2], 3.24)\n\n\n\n\n\n","category":"function"},{"location":"tutorial/varied_costs/#Working-with-varied-costs-markets","page":"Varied-costs problem","title":"Working with varied-costs markets","text":"","category":"section"},{"location":"tutorial/varied_costs/","page":"Varied-costs problem","title":"Varied-costs problem","text":"The varied-costs market is the general case of the college application problem in which schools may have different application costs. The optimization problem is as follows:","category":"page"},{"location":"tutorial/varied_costs/","page":"Varied-costs problem","title":"Varied-costs problem","text":"    beginalign*\n    textmaximizequad  v(mathcalX) = sum_jin mathcalX Bigl( f_j t_j prod_substacki in mathcalX  i  j (1 - f_i) Bigr) \n    textsubject toquad  mathcalX subseteq mathcalC sum_jin mathcalX g_j leq H\n    endalign*","category":"page"},{"location":"tutorial/varied_costs/","page":"Varied-costs problem","title":"Varied-costs problem","text":"Let's consider a random instance of this problem with m = 15 schools. ","category":"page"},{"location":"tutorial/varied_costs/","page":"Varied-costs problem","title":"Varied-costs problem","text":"using OptimalApplication\nf = rand(15)\nt = rand(60:80, 15)     # must be integers\ng = rand(5:10, 15)      # must be integers\nH = sum(g) ÷ 2\nmkt = VariedCostsMarket(f, t, g, H)","category":"page"},{"location":"tutorial/varied_costs/#Enumeration-solver","page":"Varied-costs problem","title":"Enumeration solver","text":"","category":"section"},{"location":"tutorial/varied_costs/","page":"Varied-costs problem","title":"Varied-costs problem","text":"We can solve this problem using optimalportfolio_enumerate:","category":"page"},{"location":"tutorial/varied_costs/","page":"Varied-costs problem","title":"Varied-costs problem","text":"optimalportfolio_enumerate(mkt)","category":"page"},{"location":"tutorial/varied_costs/","page":"Varied-costs problem","title":"Varied-costs problem","text":"This function returns the optimal portfolio as a vector X of school indices, and its valuation v.","category":"page"},{"location":"tutorial/varied_costs/","page":"Varied-costs problem","title":"Varied-costs problem","text":"A somewhat better kind of enumeration the branch-and-bound scheme, implemented as optimalportfolio_branchbound:","category":"page"},{"location":"tutorial/varied_costs/","page":"Varied-costs problem","title":"Varied-costs problem","text":"optimalportfolio_branchbound(mkt)","category":"page"},{"location":"tutorial/varied_costs/","page":"Varied-costs problem","title":"Varied-costs problem","text":"We can see that the results are the same (up to permutation).","category":"page"},{"location":"tutorial/varied_costs/","page":"Varied-costs problem","title":"Varied-costs problem","text":"warning: Warning\nEnumeration algorithms have exponential time complexity and should only be used in small instances with m ≤ 20. Their intended use is to test the validity of more-efficient solvers.","category":"page"},{"location":"tutorial/varied_costs/#Fast-solver","page":"Varied-costs problem","title":"Fast solver","text":"","category":"section"},{"location":"tutorial/varied_costs/","page":"Varied-costs problem","title":"Varied-costs problem","text":"OptimalApplication implements a dynamic program that solves this problem exactly in pseudopolynomial-time. The computation time of optimalportfolio_dynamicprogram depends on the values of g_j; smaller is better.","category":"page"},{"location":"tutorial/varied_costs/","page":"Varied-costs problem","title":"Varied-costs problem","text":"optimalportfolio_dynamicprogram(mkt)","category":"page"},{"location":"tutorial/varied_costs/#Approximation-algorithms","page":"Varied-costs problem","title":"Approximation algorithms","text":"","category":"section"},{"location":"tutorial/varied_costs/","page":"Varied-costs problem","title":"Varied-costs problem","text":"The college application problem admits a fully polynomial-time approximation scheme (FPTAS). Given a market and tolerance varepsilon, optimalportfolio_fptas(mkt, ε): produces a solution that is guaranteed to have an objective value of at least 1 - varepsilon times the optimum, in time O(m^3  varepsilon). ","category":"page"},{"location":"tutorial/varied_costs/","page":"Varied-costs problem","title":"Varied-costs problem","text":"optimalportfolio_fptas(mkt, 0.2)","category":"page"},{"location":"tutorial/varied_costs/","page":"Varied-costs problem","title":"Varied-costs problem","text":"warning: Warning\nSetting the tolerance ε to an extremely small value requires lots of system memory.","category":"page"},{"location":"tutorial/varied_costs/#Heuristic-algorithms","page":"Varied-costs problem","title":"Heuristic algorithms","text":"","category":"section"},{"location":"tutorial/varied_costs/","page":"Varied-costs problem","title":"Varied-costs problem","text":"OptimalApplication also includes two heuristic algorithms. Although they offer no accuracy guarantee, their computation time is low and they can be quite effective in practice.","category":"page"},{"location":"tutorial/varied_costs/","page":"Varied-costs problem","title":"Varied-costs problem","text":"optimalportfolio_greedy adds schools in descending order by f_j t_j / g_j until the budget is exhausted.","category":"page"},{"location":"tutorial/varied_costs/","page":"Varied-costs problem","title":"Varied-costs problem","text":"optimalportfolio_greedy(mkt)","category":"page"},{"location":"tutorial/varied_costs/","page":"Varied-costs problem","title":"Varied-costs problem","text":"optimalportfolio_simulatedannealing starts with the greedy solution, then searches locally according to a randomized scheme for an improvement.","category":"page"},{"location":"tutorial/varied_costs/","page":"Varied-costs problem","title":"Varied-costs problem","text":"optimalportfolio_simulatedannealing(mkt)","category":"page"},{"location":"public_api/market_io/#Market-input-and-output","page":"Market I/O","title":"Market input and output","text":"","category":"section"},{"location":"public_api/market_io/","page":"Market I/O","title":"Market I/O","text":"Market\nSameCostsMarket\nVariedCostsMarket\nvaluation","category":"page"},{"location":"public_api/market_io/#OptimalApplication.Market","page":"Market I/O","title":"OptimalApplication.Market","text":"Contains information about a college application market.\n\n\n\n\n\n","category":"type"},{"location":"public_api/market_io/#OptimalApplication.SameCostsMarket","page":"Market I/O","title":"OptimalApplication.SameCostsMarket","text":"SameCostsMarket{U<:Real}\n\nContains information about a college application market with identical application costs.\n\nConstructors\n\nSameCostsMarket(f::Vector{Float64}, t::Vector{U}, h::Integer) where {U<:Real}\n\nConstruct the SameCostsMarket defined by admissions probabilities f, utility values t, and application limit h.\n\njulia> SameCostsMarket([0.3, 0.2, 0.05], [2, 3, 4], 2)\nSameCostsMarket{Int64}(3, [0.3, 0.2, 0.05], [2, 3, 4], 2, [0.6, 0.6000000000000001, 0.2], [0.7, 0.8, 0.95], [1, 2, 3])\n\nSameCostsMarket(m)\n\nGenerate a random SameCostsMarket{Int} with m schools. f and t correlate negatively, mimicking a realistic market.\n\nInternal API\n\nThis type contains the following fields:\n\nm: Number of schools\nf: Vector of admissions probabilities\nt: Vector of utility values (must be sorted)\nh: Number of schools student is allowed to apply to\nft = f .* t\nomf = 1 .- f\nperm: How the input data were permuted to sort by t\n\nU is the eltype of t.\n\nInternally, the schools must be indexed so that t is sorted ascending. perm allows the original order to be recovered:\n\njulia> f = rand(3); t = rand(3); mkt = SameCostsMarket(f, t, 2);\n\njulia> t[mkt.perm] == mkt.t\ntrue\n\njulia> t == mkt.t[invperm(mkt.perm)]\ntrue\n\n\n\n\n\n","category":"type"},{"location":"public_api/market_io/#OptimalApplication.VariedCostsMarket","page":"Market I/O","title":"OptimalApplication.VariedCostsMarket","text":"VariedCostsMarket\n\nContains information about a college application market with varied application costs.\n\nConstructors\n\nVariedCostsMarket(f::Vector{Float64}, t::Vector{Int}, g::Vector{Int}, H::Int)\n\nConstruct the VariedCostsMarket defined by admissions probabilities f, utility values t, application costs g, and application limit h.\n\njulia> mkt = VariedCostsMarket([0.1, 0.5, 0.3, 0.1], [12, 3, 4, 13], [1, 2, 1, 1], 3)\nVariedCostsMarket(4, [0.5, 0.3, 0.1, 0.1], [3, 4, 12, 13], [2, 1, 1, 1], 3, [1.5, 1.2, 1.2000000000000002, 1.3], [0.5, 0.7, 0.9, 0.9], [2, 3, 1, 4])\n\nVariedCostsMarket(m)\n\nGenerate a random SameCostsMarket{Int} with m schools. f and t correlate negatively, mimicking a realistic market. Entries of g are random integers in 5:10.\n\nInternal API\n\nThis type contains the following fields:\n\nm: Number of schools\nf: Vector of admissions probabilities\nt: Vector of utility values (must be sorted)\ng: Vector of application costs\nH: Budget to spend on applications\nft = f .* t\nomf = 1 .- f\nperm: How the input data were permuted to sort by t\n\nInternally, the schools must be indexed so that t is sorted ascending. perm allows the original order to be recovered:\n\njulia> f = rand(3); t = rand(3); g = rand(5:10, 3); mkt = VariedCostsMarket(f, t, g, 10);\n\njulia> t[mkt.perm] == mkt.t\ntrue\n\njulia> t == mkt.t[invperm(mkt.perm)]\ntrue\n\n\n\n\n\n","category":"type"},{"location":"public_api/market_io/#OptimalApplication.valuation","page":"Market I/O","title":"OptimalApplication.valuation","text":"valuation(X, mkt)\n\nReturn the valuation of the portfolio X for the market mkt, which may be either  a SameCostsMarket or a VariedCostsMarket.\n\nSameCostsMarket example (h = 3 is irrelevant):\n\njulia> mkt = SameCostsMarket([0.1, 0.5, 0.3, 0.1], [12, 3, 4, 13], 3);\n\njulia> round(valuation([1, 4], mkt), digits=2) # expected utility when applying to schools 1 and 4\n2.38\n\nVariedCostsMarket example (g = [1, 2, 1, 1] and H = 4 are irrelevant):\n\njulia> mkt = VariedCostsMarket([0.1, 0.5, 0.3, 0.1], [12, 3, 4, 13], [1, 2, 1, 1], 4);\n\njulia> round(valuation([1, 4], mkt), digits=2) # expected utility when applying to schools 1 and 4\n2.38\n\n\n\n\n\n","category":"function"},{"location":"tutorial/market_io/#Market-input-and-output","page":"Market I/O","title":"Market input and output","text":"","category":"section"},{"location":"tutorial/market_io/","page":"Market I/O","title":"Market I/O","text":"To read in an instance of the same-costs problem, use SameCostsMarket(f, t, h):","category":"page"},{"location":"tutorial/market_io/","page":"Market I/O","title":"Market I/O","text":"using OptimalApplication\nmkt = SameCostsMarket(\n    [0.5, 0.1, 0.9, 0.7],   # admissions probabilities\n    [12, 20, 1, 3],         # utility values\n    2                       # limit on number of schools to apply to\n)","category":"page"},{"location":"tutorial/market_io/","page":"Market I/O","title":"Market I/O","text":"The output is an instance of SameCostsMarket{U}, which contains the input data along with certain preliminary calculations, such as 1 .- f, that are used repeatedly in the solution algorithms. Here  U == eltype(t) and we require U <: Real.","category":"page"},{"location":"tutorial/market_io/","page":"Market I/O","title":"Market I/O","text":"note: Note\nTo improve the performance of the solvers, the schools are sorted in ascending order by t when the market is constructed. Therefore, mkt = SameCostsMarket(f, t, h); mkt.t == t is not necessarily true. The input data can be recovered using mkt.t == t[mkt.perm]; see  SameCostsMarket(f, t, h).","category":"page"},{"location":"tutorial/market_io/","page":"Market I/O","title":"Market I/O","text":"To read in an instance of the varied-costs problem, use VariedCostsMarket(f, t, g, H):","category":"page"},{"location":"tutorial/market_io/","page":"Market I/O","title":"Market I/O","text":"mkt = VariedCostsMarket(\n    [0.5, 0.1, 0.9, 0.7],   # admissions probabilities\n    [12, 20, 1, 3],         # utility values\n    [3, 4, 7, 2],           # application fees\n    8                       # budget to spend on applications\n)","category":"page"},{"location":"tutorial/market_io/","page":"Market I/O","title":"Market I/O","text":"The output is an instance of VariedCostsMarket. Both SameCostsMarket and VariedCostsMarket are subtypes  of the abstract type Market.","category":"page"},{"location":"tutorial/market_io/","page":"Market I/O","title":"Market I/O","text":"note: Note\nVariedCostsMarket supports only t with integer eltype. Because the objective-function is linear in t, to work with float data, first multiply by the least common denominator. ","category":"page"},{"location":"tutorial/market_io/","page":"Market I/O","title":"Market I/O","text":"If the application fees share a common divisor, the solvers will perform more effectively if you divide by it. For example, fees of [90, 80, 90, 70] with budget 160 is equivalent to fees of [9, 8, 9, 7] with budget 16. ","category":"page"},{"location":"tutorial/market_io/","page":"Market I/O","title":"Market I/O","text":"OptimalApplication also provides the convenience functions Market(f, t, h) and Market(f, t, g, H)](@ref), which behave just like those above.","category":"page"},{"location":"public_api/enumeration_solvers/#Enumeration-solvers","page":"Enumeration solvers","title":"Enumeration solvers","text":"","category":"section"},{"location":"public_api/enumeration_solvers/","page":"Enumeration solvers","title":"Enumeration solvers","text":"optimalportfolio_enumerate\noptimalportfolio_branchbound","category":"page"},{"location":"public_api/enumeration_solvers/#OptimalApplication.optimalportfolio_enumerate","page":"Enumeration solvers","title":"OptimalApplication.optimalportfolio_enumerate","text":"optimalportfolio_enumerate(mkt::SameCostsMarket) -> X, v\n\nProduce the optimal portfolio for SameCostsMarket defined by mkt. Very slow; you probably should use applicationorder_list or applicationorder_heap instead.\n\njulia> mkt = SameCostsMarket([0.2, 0.5, 0.1, 0.6, 0.1], [1, 4, 9, 1, 8], 2);\n\njulia> optimalportfolio_enumerate(mkt)\n([2, 3], 2.7)\n\n\n\n\n\noptimalportfolio_enumerate(mkt::VariedCostsMarket) -> X, v\n\nProduce the optimal portfolio for VariedCostsMarket defined by mkt. Very slow; you probably should use optimalportfolio_dynamicprogram or optimalportfolio_branchbound or optimalportfolio_simulatedannealing instead.\n\njulia> mkt = VariedCostsMarket([0.2, 0.5, 0.1, 0.6, 0.1], [1, 4, 9, 1, 8], [2, 4, 2, 5, 1], 8);\n\njulia> optimalportfolio_enumerate(mkt)\n([2, 5, 3], 3.24)\n\n\n\n\n\n","category":"function"},{"location":"public_api/enumeration_solvers/#OptimalApplication.optimalportfolio_branchbound","page":"Enumeration solvers","title":"OptimalApplication.optimalportfolio_branchbound","text":"optimalportfolio_branchbound(mkt::VariedCostsMarket; maxit=10000, verbose=false) -> X, v\n\nUse the branch-and-bound algorithm to produce the optimal portfolio for the VariedCostsMarket defined by mkt. Intractable for large markets. \n\njulia> mkt = VariedCostsMarket([0.2, 0.5, 0.1, 0.6, 0.1], [1, 4, 9, 1, 8], [2, 4, 2, 5, 1], 8);\n\njulia> optimalportfolio_branchbound(mkt)\n([2, 3, 5], 3.24)\n\n\n\n\n\n","category":"function"},{"location":"tutorial/same_costs/#Working-with-same-costs-markets","page":"Same-costs problem","title":"Working with same-costs markets","text":"","category":"section"},{"location":"tutorial/same_costs/","page":"Same-costs problem","title":"Same-costs problem","text":"The same-costs market is a special case of the college application problem in which all schools have the same application cost. In this case, the budget constraint reduces to a limit on the number of schools to which you can apply:","category":"page"},{"location":"tutorial/same_costs/","page":"Same-costs problem","title":"Same-costs problem","text":"    beginalign*\n    textmaximizequad  v(mathcalX) = sum_jin mathcalX Bigl( f_j t_j prod_substacki in mathcalX  i  j (1 - f_i) Bigr) \n    textsubject toquad  mathcalX subseteq mathcalC lvert mathcalX rvert leq h\n    endalign*","category":"page"},{"location":"tutorial/same_costs/","page":"Same-costs problem","title":"Same-costs problem","text":"Let's consider a random instance of this problem with m = 15 schools and h = 5. ","category":"page"},{"location":"tutorial/same_costs/","page":"Same-costs problem","title":"Same-costs problem","text":"using OptimalApplication\nf = rand(15)\nt = rand(15)\nmkt = SameCostsMarket(f, t, 5)","category":"page"},{"location":"tutorial/same_costs/#Enumeration-solver","page":"Same-costs problem","title":"Enumeration solver","text":"","category":"section"},{"location":"tutorial/same_costs/","page":"Same-costs problem","title":"Same-costs problem","text":"We can solve this problem using optimalportfolio_enumerate:","category":"page"},{"location":"tutorial/same_costs/","page":"Same-costs problem","title":"Same-costs problem","text":"X, v = optimalportfolio_enumerate(mkt)","category":"page"},{"location":"tutorial/same_costs/","page":"Same-costs problem","title":"Same-costs problem","text":"This function returns the optimal portfolio as a vector X of school indices, and its valuation v.","category":"page"},{"location":"tutorial/same_costs/","page":"Same-costs problem","title":"Same-costs problem","text":"warning: Warning\nEnumeration algorithms have exponential time complexity and should only be used in small instances with m ≤ 20. Their intended use is to test the validity of more-efficient solvers.","category":"page"},{"location":"tutorial/same_costs/#nestednessproperty","page":"Same-costs problem","title":"The nestedness property","text":"","category":"section"},{"location":"tutorial/same_costs/","page":"Same-costs problem","title":"Same-costs problem","text":"The optimal solutions for this special case of the college application problem satisfy a property called the nestedness property. When mathcalX_h denotes the optimal solution for a given market when the limit is h, we have","category":"page"},{"location":"tutorial/same_costs/","page":"Same-costs problem","title":"Same-costs problem","text":"mathcalX_1 subset mathcalX_2  subset dots subset mathcalX_m","category":"page"},{"location":"tutorial/same_costs/","page":"Same-costs problem","title":"Same-costs problem","text":"This means that there exists a permutation of the schools X such that X[1:h] gives the optimal portfolio of size h for all h.","category":"page"},{"location":"tutorial/same_costs/#Efficient-solvers","page":"Same-costs problem","title":"Efficient solvers","text":"","category":"section"},{"location":"tutorial/same_costs/","page":"Same-costs problem","title":"Same-costs problem","text":"The functions applicationorder_list and applicationorder_heap produces this permutation X as well as a vector V indicating the utility associated with each optimum. When mkt.h < mkt.m, only the first mkt.h entries are computed. To verify the correctness of the nestedness property, however, let's make a copy of our market with mkt.h == 15 and check that the first h = 5 entries agree (up to permutation).","category":"page"},{"location":"tutorial/same_costs/","page":"Same-costs problem","title":"Same-costs problem","text":"mkt_long = SameCostsMarket(f, t, 15)\nX_long, V_long = applicationorder_list(mkt)\nX_long[1:5], X","category":"page"},{"location":"tutorial/same_costs/","page":"Same-costs problem","title":"Same-costs problem","text":"The valuation should also be the same:","category":"page"},{"location":"tutorial/same_costs/","page":"Same-costs problem","title":"Same-costs problem","text":"V_long[5], v","category":"page"},{"location":"tutorial/same_costs/","page":"Same-costs problem","title":"Same-costs problem","text":"The output of applicationorder_heap is identical to that of applicationorder_list but the function uses a different internal data structure. In the vast majority of cases, the list version is faster.","category":"page"},{"location":"examples/valuation/#Planets-example","page":"Valuation curve","title":"Planets example","text":"","category":"section"},{"location":"examples/valuation/","page":"Valuation curve","title":"Valuation curve","text":"Let's consider a same-costs market and plot the valuation curve, which gives v(mathcalX_h) as a function of h, where mathcalX_h is the optimal portfolio for each h between 0 and m.","category":"page"},{"location":"examples/valuation/","page":"Valuation curve","title":"Valuation curve","text":"The market data is as follows:","category":"page"},{"location":"examples/valuation/","page":"Valuation curve","title":"Valuation curve","text":"beginalign\nf = (039 033 024 024 005 003 01 012) \nt = (200 250 300 350 400 450 500 550)\nendalign","category":"page"},{"location":"examples/valuation/","page":"Valuation curve","title":"Valuation curve","text":"Let's import the PrettyTables and UnicodePlots libraries for visualizations and read in the market data.","category":"page"},{"location":"examples/valuation/","page":"Valuation curve","title":"Valuation curve","text":"using PrettyTables\nusing UnicodePlots\nusing OptimalApplication\nf = [0.39, 0.33, 0.24, 0.24, 0.05, 0.03, 0.1, 0.12]\nt = [200, 250, 300, 350, 400, 450, 500, 550]\nmkt = Market(f, t, 8)","category":"page"},{"location":"examples/valuation/","page":"Valuation curve","title":"Valuation curve","text":"Notice how we set h = m = 8 because we know from the nestedness property that by calling applicationorder_list(mkt), we can get a permutation of the schools that encodes all the optimal portfolios. Let's do that now:","category":"page"},{"location":"examples/valuation/","page":"Valuation curve","title":"Valuation curve","text":"X, V = applicationorder_list(mkt, verbose=true)","category":"page"},{"location":"examples/valuation/","page":"Valuation curve","title":"Valuation curve","text":"mathcalX_h is given by the first h entries of X, and v(mathcalX_h) is V[h]. Let's plot V.","category":"page"},{"location":"examples/valuation/","page":"Valuation curve","title":"Valuation curve","text":"lineplot(0:mkt.m, vcat(0, V), xlabel=\"h\", ylabel=\"v\")","category":"page"},{"location":"examples/valuation/","page":"Valuation curve","title":"Valuation curve","text":"This graph is called a valuation curve. It gives the applicant's utility as a function of her budget h. The curve's concave shape is guaranteed by the nestedness property. ","category":"page"},{"location":"examples/valuation/","page":"Valuation curve","title":"Valuation curve","text":"An alternative way of looking at the optimal portfolios is to take the inverse-permutation of X. The value of invperm(X)[j] gives the smallest value of h for which school j is in mathcalX_h. This is called the schools priority number; if a school has priority number 1, it is the first school you should apply to, and so on. ","category":"page"},{"location":"examples/valuation/","page":"Valuation curve","title":"Valuation curve","text":"This format lets us create the following table:","category":"page"},{"location":"examples/valuation/","page":"Valuation curve","title":"Valuation curve","text":"priority = invperm(X)\npretty_table(\n    Any[f t priority V[priority]], \n    header = [\"f\", \"t\", \"priority\", \"valuation\"],\n    header_crayon = crayon\"bold yellow\",\n)","category":"page"},{"location":"#OptimalApplication","page":"Home","title":"OptimalApplication","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"... is a Julia package for solving the college application problem, which is the following optimization problem:","category":"page"},{"location":"","page":"Home","title":"Home","text":"    beginalign*\n    textmaximizequad  v(mathcalX) =  operatornameEBiglmaxbigr0\n    maxt_j Z_j  j in mathcalXbigrBigr \n    textsubject toquad  mathcalX subseteq mathcalC sum_jin mathcalX g_j leq H\n    endalign*","category":"page"},{"location":"","page":"Home","title":"Home","text":"Here mathcalC =  1 dots m is a collection of schools, t_j is the utility associated with attending school j, and Z_j is a random Bernoulli variable with probability f_j. The problem is to decide on a subset of schools mathcalX to apply to, given g_j, the application fee for school j, and H the total budget to spend on applications. v(mathcalX) represents the expected utility associated with the application portfolio mathcalX, and by indexing the schools such that each t_j leq t_j+1, it can be written in the following closed form: ","category":"page"},{"location":"","page":"Home","title":"Home","text":"    v(mathcalX) = sum_jin mathcalX Bigl( f_j t_j prod_substacki in mathcalX  i  j (1 - f_i) Bigr)","category":"page"},{"location":"","page":"Home","title":"Home","text":"The general college application problem is referred to in this package as the varied-costs college application problem. A special case, called the same-costs problem, occurs when each g_j = 1. Then the budget constraint becomes a simple cardinality constraint lvert mathcalX rvert leq h, where h is written in lowercase to reflect the change in  meaning. ","category":"page"},{"location":"","page":"Home","title":"Home","text":"The same-costs college application problem is solvable in polynomial time by a greedy algorithm, whereas the varied-costs problem is NP-complete (as can be shown by a reduction from the knapsack problem).","category":"page"},{"location":"","page":"Home","title":"Home","text":"This package provides:","category":"page"},{"location":"","page":"Home","title":"Home","text":"Types SameCostsMarket and VariedCostsMarket for holding problem instance data.\nVarious exact algorithms for SameCostsMarkets\nExact algorithms for VariedCostsMarkets\nAn approximation algorithm for VariedCostsMarkets\nHeuristic algorithms for VariedCostsMarkets","category":"page"},{"location":"","page":"Home","title":"Home","text":"Pages = [\n    \"man/same_costs.md\",\n    \"man/varied_costs.md\",\n]\nDepth = 1","category":"page"}]
}
