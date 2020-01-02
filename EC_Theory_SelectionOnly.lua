dofile("Statistics.lua")
dofile("Bitmap.lua")
dofile("Graphics.lua")
dofile("EC_Common.lua")
dofile("GA_Common.lua")

local RUNS                          = 400
local MAX_POP_SIZE                  = 200
local MAX_TRUNCATION_FACTOR         = 25

local GENOME_BITS                   = 32
local GENOME_MAX_VALUE              = math.pow(2, GENOME_BITS) - 1
local MAX_FITNESS                   = 100
local FITNESS_FUNCTION              = false

local IMAGE_WIDTH                   = 1000
local IMAGE_HEIGHT                  = 1000
local IMAGE_FILENAME_FIGURE_6_1     = "EC/Figure6_1.bmp"
local IMAGE_FILENAME_FIGURE_6_2     = "EC/Figure6_2.bmp"
local IMAGE_FILENAME_FIGURE_6_3     = "EC/Figure6_3.bmp"
local IMAGE_FILENAME_FIGURE_6_4     = "EC/Figure6_4.bmp"
local IMAGE_FILENAME_FIGURE_6_5     = "EC/Figure6_5.bmp"
local IMAGE_FILENAME_FIGURE_6_6     = "EC/Figure6_6.bmp"
local IMAGE_FILENAME_FIGURE_6_7     = "EC/Figure6_7.bmp"
local IMAGE_FILENAME_FIGURE_6_8     = "EC/Figure6_8.bmp"
local IMAGE_FILENAME_FIGURE_6_9     = "EC/Figure6_9.bmp"
local IMAGE_FILENAME_FIGURE_6_10    = "EC/Figure6_10.bmp"
local IMAGE_FILENAME_FIGURE_6_11    = "EC/Figure6_11.bmp"
local IMAGE_FILENAME_FIGURE_6_12    = "EC/Figure6_12.bmp"
local IMAGE_FILENAME_FIGURE_6_13    = "EC/Figure6_13.bmp"
local IMAGE_FILENAME_FIGURE_6_14    = "EC/Figure6_14.bmp"
local IMAGE_FILENAME_FIGURE_6_15    = "EC/Figure6_15.bmp"
local IMAGE_FILENAME_FIGURE_6_16    = "EC/Figure6_16.bmp"
local IMAGE_FILENAME_FIGURE_6_17    = "EC/Figure6_17.bmp"
local IMAGE_FILENAME_FIGURE_6_18    = "EC/Figure6_18.bmp"

local s_AssignedFitness = {}

local function AssignUniqueFitness(pop)
  local used, used_count = {}, 0
  for i, ind in ipairs(pop) do
    ind.fitness = math.random(0, MAX_FITNESS - 1)
    while used[ind.fitness] do
      ind.fitness = math.random(0, MAX_FITNESS - 1)
    end
    used[ind.fitness] = true
    used_count = used_count + 1
    if used_count == MAX_FITNESS then
      used, used_count = {}, 0
    end
    s_AssignedFitness[ind.chrom] = ind.fitness
  end
  local best = pop[math.random(1, #pop)]
  best.fitness = MAX_FITNESS
  s_AssignedFitness[best.chrom] = best.fitness
end

local function GetAssignedFitness(chrom)
  return s_AssignedFitness[chrom]
end

local function EvalPop(pop)
  local min_fitness, max_fitness
  local total_fitness = 0.0
  local best_ind
  for _, ind in ipairs(pop) do
    local fitness = FITNESS_FUNCTION(ind.chrom)
    ind.fitness = fitness
    min_fitness = (not min_fitness or fitness < min_fitness) and fitness or min_fitness
    max_fitness = (not max_fitness or fitness > max_fitness) and fitness or max_fitness
    total_fitness = total_fitness + fitness
    best_ind = (not best_ind or fitness > best_ind.fitness) and ind or best_ind
  end
  pop.min_fitness, pop.max_fitness, pop.avg_fitness = min_fitness, max_fitness, total_fitness / #pop
  pop.spread_fitness = max_fitness - min_fitness
  pop.best_ind = best_ind
end

local function GenUniqueInitPop(size, on_birth)
  local pop = {}
  local used = {}
  while #pop < size do
    local chrom
    while not chrom or used[chrom] do
      chrom = math.random(0, GENOME_MAX_VALUE)
    end
    used[chrom] = true
    table.insert(pop, {chrom = chrom, birth = #pop + 1})
    if on_birth then
      EvalPop(pop)
      on_birth(pop, #pop)
    end
  end

  return pop
end

-- TODO: this should catch the case with >1 individual left too, e.g. last 10 gens without change - return the gen which started the "no-change"
local function IsPopConverged(pop)
  local unique = pop[1].chrom
  for _, ind in ipairs(pop) do
    if ind.chrom ~= unique then
      return false
    end
  end
  
  return true
end

local function UniformStochasticSelection(parents, size)
  local offsprings = {}
  for i = 1, size do
    local parent = parents[math.random(1, #parents)]
    table.insert(offsprings, parent)
  end
  
  return offsprings
end

local function GetMergedSorted(parents, offsprings)
  local sorted = {}
  for _, ind in ipairs(parents) do
    table.insert(sorted, ind)
  end
  for _, ind in ipairs(offsprings) do
    table.insert(sorted, ind)
  end
  table.sort(sorted, function(a, b) return a.fitness > b.fitness end)
  
  return sorted
end

local function TruncateSelection(parents, size, offsprings)
  offsprings = offsprings or {}
  
  local pool = GetMergedSorted(parents, offsprings)
  offsprings = {}
  while #offsprings < size do
    local parent = pool[math.random(1, size)]
    table.insert(offsprings, parent)
  end
  
  return offsprings
end

local function GetLinearRankingProb(parents)
  local m = #parents
  local epsilon = 1.0 / (m * m)
  table.sort(parents, function(a, b) return a.fitness > b.fitness end)
  
  local prob, sum = {}, 0.0
  for i = 1, m do
    local p = (2.0 / m - epsilon) - (2.0 / m - 2 * epsilon) * (i - 1) / (m - 1)
    sum = sum + p
    prob[i] = sum
  end
  assert(math.abs(sum - 1.0) < 0.0001)
  
  return prob
end

local function LinearRankingSelection(parents, size)
  local prob = GetLinearRankingProb(parents)
  local offsprings = {}
  while #offsprings < size do
    local rand = math.random()
    local i = table.binary_search(prob, rand)
    table.insert(offsprings, parents[i])
  end
  
  return offsprings
end

local function BinaryTournamentSelection(parents, size)
  local m = #parents
  local offsprings = {}
  while #offsprings < size do
    local ind1 = parents[math.random(1, m)]
    local ind2 = parents[math.random(1, m)]
    table.insert(offsprings, (ind1.fitness > ind2.fitness) and ind1 or ind2)
  end
  
  return offsprings
end

local function DeterministicSelection(parents, size)
  local brood_size = size // #parents
  local offsprings = {}
  for _, parent in ipairs(parents) do
    for i = 1, brood_size do
      table.insert(offsprings, parent)
    end
  end
  while #offsprings < size do
    table.insert(offsprings, parents[size - #offsprings])
  end
  
  return offsprings
end

local function GetFitnessProportionalProb(parents)
  local prob, total_prob = {}, 0.0
  for i, ind in ipairs(parents) do
    total_prob = total_prob + ind.fitness
    prob[i] = total_prob
  end
  prob.total = total_prob
  
  return prob
end

local function FitnessProportionalSelection(parents, size)
  local prob = GetFitnessProportionalProb(parents)
  local offsprings = {}
  while #offsprings < size do
    local rand = math.random() * prob.total
    local i = table.binary_search(prob, rand)
    table.insert(offsprings, parents[i])
  end
  
  return offsprings
end

local function SUS_WheelOfFortune(parents, size, prob, slot_size)
  local offsprings = {}
  local rand = math.random() * slot_size
  local i = table.binary_search(prob, rand)
  while #offsprings < size do
    table.insert(offsprings, parents[i])
    rand = rand + slot_size
    while i < #parents and prob[i] < rand do
      i = i + 1
    end
  end
  
  return offsprings
end

local function SUS_FitnessProportionalSelection(parents, size)
  local prob = GetFitnessProportionalProb(parents)
  local slot_size = prob.total / size
  
  return SUS_WheelOfFortune(parents, size, prob, slot_size)  
end

local function SUS_LinearRanking(parents, size)
  local prob = GetLinearRankingProb(parents)
  local slot_size = 1.0 / size
  
  return SUS_WheelOfFortune(parents, size, prob, slot_size)
end

local function EV(m, n, selection_parent, selection_survival, overlapping)
  local parents = GenUniqueInitPop(m)
  AssignUniqueFitness(parents)
  local gen = 0
  while not IsPopConverged(parents) do
    local offsprings = selection_parent(parents, n)
    if m == n and not selection_survival then
      parents = offsprings
    else
      if overlapping then
        for _, parent in ipairs(parents) do
          table.insert(offsprings, parent)
        end
      end
      parents = (selection_survival or UniformStochasticSelection)(offsprings, m)
    end
    --EvalPop(parents)
    gen = gen + 1
  end
  
  return gen, parents
end

-- TODO: this can be implemented via EV()
local function GA(m, selection_parent)
  local parents = GenUniqueInitPop(m)
  AssignUniqueFitness(parents)
  local gen = 0
  while not IsPopConverged(parents) do
    parents = selection_parent(parents, m)
    gen = gen + 1
  end
  
  return gen, parents
end

local s_Seeds = {}
for run = 1, MAX_POP_SIZE do
  s_Seeds[run] = math.random(1, 100000)
end

local function CompareNaturalSelection(filename, filename2)
  local graphs = {name_x = "Population Size m", name_y = "Average Generations to Fixed-Point Convergence", funcs = {}}
  local name = string.format("fixed(m)")
  local points = {color = RGB_GREEN, stddev_start = 1, stddev_interval = 50, stddev_color = RGB_CYAN}
  graphs.funcs[name] = points
  
  local graphs2 = {name_x = "Population Size m", name_y = "Average Fixed-Point Fitness", funcs = {}}
  local points2 = {color = RGB_GREEN}
  graphs2.funcs[name] = points2

  for m = 10, MAX_POP_SIZE do
    local total_gens, total_fitness = 0, 0
    local values = {}
    for run = 1, RUNS do
      local gens, pop = EV(m, m, UniformStochasticSelection)
      total_gens = total_gens + gens
      total_fitness = total_fitness + pop[1].fitness
      table.insert(values, gens)
    end
    table.insert(points, {x = m, y = total_gens, values = ((#points + 1) - points.stddev_start) % points.stddev_interval == 0 and values})
    table.insert(points2, {x = m, y = total_fitness})
    print(string.format("Population Size: %d, Average generations/fitness to fixed point: %.2f/%.2f", m, total_gens / RUNS, total_fitness / RUNS))
  end
  
  NormalizeGraphs(graphs, RUNS)
  NormalizeGraphs(graphs2, RUNS)
  
  local bmp = Bitmap.new(IMAGE_WIDTH, IMAGE_HEIGHT, RGB_BLACK)
  DrawGraphs(bmp, graphs, {int_x = true, skip_KP = true, min_x = 0, min_y = 0})
  bmp:WriteBMP(filename)
  local bmp = Bitmap.new(IMAGE_WIDTH, IMAGE_HEIGHT, RGB_BLACK)
  DrawGraphs(bmp, graphs2, {int_x = true, skip_KP = true, min_x = 0, min_y = 0})
  bmp:WriteBMP(filename2)
end

local function CompareTruncFactor(m, filename, filename2)
  local graphs = {name_x = "Truncation Factor k", name_y = "Average Generations to Fixed-Point Convergence", funcs = {}}
  local name_infinite = string.format("Infinite Population Model for m=%d", m)
  local name_finite = string.format("Finite Stochastic Population Model for m=%d", m)
  local points_infinite = {color = RGB_GREEN}
  local points_finite = {color = RGB_CYAN}
  graphs.funcs[name_infinite] = points_infinite
  graphs.funcs[name_finite] = points_finite
  
  local graphs2 = {name_x = "Truncation Factor k", name_y = "Average Fixed-Point Fitness", funcs = {}}
  local name2_infinite = "Infinite-Population homogenous fixed-point"
  local name2_finite = "Finite-Population homogenoues fixed-point"
  local points2_infinite = {color = RGB_GREEN}
  local points2_finite = {color = RGB_CYAN}
  graphs2.funcs[name2_infinite] = points2_infinite
  graphs2.funcs[name2_finite] = points2_finite
  
  for k = 0, MAX_TRUNCATION_FACTOR do
    local gens_infinite = math.ceil(math.log(m) / math.log(m / k))
    local gens_finite, fitness_finite
    if k == 0 then
      gens_finite, fitness_finite = 0, MAX_FITNESS
    else
      math.randomseed(s_Seeds[m])
      gens_finite, fitness_finite = 0, 0
      for run = 1, RUNS do
        local gens, pop = EV(m, k, TruncateSelection)
        gens_finite = gens_finite + gens
        fitness_finite = fitness_finite + pop[1].fitness
      end
      gens_finite, fitness_finite = gens_finite / RUNS, fitness_finite / RUNS
    end
    table.insert(points_infinite, {x = k, y = gens_infinite})
    table.insert(points_finite, {x = k, y = gens_finite})
    table.insert(points2_infinite, {x = k, y = MAX_FITNESS})
    table.insert(points2_finite, {x = k, y = fitness_finite})
    print(k, gens_finite, fitness_finite)
  end
  
  local bmp = Bitmap.new(IMAGE_WIDTH, IMAGE_HEIGHT, RGB_BLACK)
  DrawGraphs(bmp, graphs, {int_x = true, skip_KP = true, min_x = 0, min_y = 0, div_x = 5, div_y = 5, max_y = 10})
  bmp:WriteBMP(filename)
  local bmp = Bitmap.new(IMAGE_WIDTH, IMAGE_HEIGHT, RGB_BLACK)
  DrawGraphs(bmp, graphs2, {int_x = true, skip_KP = true, min_x = 0, min_y = 96, div_x = 5, div_y = 5})
  bmp:WriteBMP(filename2)
end

local function TakeOverCurves(m, filename, filename2, filename3, filename4)
  local graphs = {name_x = "Generation", name_y = "TakeOver Proportion", funcs = {}}
  local name_linear = string.format("Linear Ranking with m=%d and epsilon=1/(m*m)", m)
  local points_linear = {color = RGB_GREEN}
  local points_trunc1 = {color = RGB_CYAN}
  local points_trunc5 = {color = RGB_WHITE}
  local name_trunc1 = string.format("Truncation Selection for m=%d and k=1", m)
  local name_trunc5 = string.format("Truncation Selection for m=%d and k=5", m)
  graphs.funcs[name_linear] = points_linear
  graphs.funcs[name_trunc1] = points_trunc1
  graphs.funcs[name_trunc5] = points_trunc5
  
  local epsilon = 1.0 / (m * m)
  local Pt_linear, Pt_trunc1, Pt_trunc5 = 1.0 / m, 1.0 / m, 1.0 / m
  table.insert(points_linear, {x = 0, y = Pt_linear})
  for gen = 1, 12 do
    Pt_linear = 2 * Pt_linear - Pt_linear * Pt_linear
    table.insert(points_linear, {x = gen, y = Pt_linear})
  end
  table.insert(points_trunc1, {x = 0, y = Pt_trunc1})
  table.insert(points_trunc5, {x = 0, y = Pt_trunc5})
  for gen = 0.1, 12, 0.1 do
    Pt_trunc1 = clamp((1.0 / m) * math.pow(m / 1, gen), 0, 1)
    Pt_trunc5 = clamp((1.0 / m) * math.pow(m / 5, gen), 0, 1)
    table.insert(points_trunc1, {x = gen, y = Pt_trunc1})
    table.insert(points_trunc5, {x = gen, y = Pt_trunc5})
  end
  
  local bmp = Bitmap.new(IMAGE_WIDTH, IMAGE_HEIGHT, RGB_BLACK)
  DrawGraphs(bmp, graphs, {int_x = true, skip_KP = true, min_x = 0, min_y = 0, div_x = 6, div_y = 5, max_y = 1.0})
  bmp:WriteBMP(filename)
  
  local graphs2 = {name_x = "Population Size m", name_y = "Average Generations to Fixed-Point Convergence", funcs = {}}
  local name_infinite = "Infinite-Population Model"
  local name_finite = "Finite-Population Model"
  local points_infinite = {color = RGB_GREEN}
  local points_finite = {color = RGB_CYAN}
  graphs2.funcs[name_infinite] = points_infinite
  graphs2.funcs[name_finite] = points_finite

  local graphs3 = {name_x = "Parent=Offspring Population Size", name_y = "Average Generations to Fixed-Point Convergence", funcs = {}}
  local name_binary = "Observed Binary Tournament"
  local name_ranking = "Observed Linear Ranking"
  local points_binary = {color = RGB_GREEN}
  local points_ranking = {color = RGB_CYAN}
  graphs3.funcs[name_binary] = points_binary
  graphs3.funcs[name_ranking] = points_ranking
  
  local graphs4 = {name_x = "Parent=Offspring Population Size", name_y = "Average Population Fitness at Fixed-Point Convergence", funcs = {}}
  local name_binary_fitness = "Binary Tournament"
  local name_ranking_fitness = "Linear Ranking"
  local points_binary_fitness = {color = RGB_GREEN, stddev_start = 1, stddev_interval = 25}
  local points_ranking_fitness = {color = RGB_CYAN, stddev_start = 2, stddev_interval = 25}
  graphs4.funcs[name_binary_fitness] = points_binary_fitness
  graphs4.funcs[name_ranking_fitness] = points_ranking_fitness
  
  for size = 10, MAX_POP_SIZE do
    local gens_infinite = math.log(size) + math.log(math.log(size))
    table.insert(points_infinite, {x = size, y = gens_infinite})
    math.randomseed(s_Seeds[size])
    local gens_ranking, fitness_ranking, values_ranking = 0, 0.0, {}
    for i = 1, RUNS do
      local gens, pop = EV(size, size, LinearRankingSelection)
      gens_ranking = gens_ranking + gens
      fitness_ranking = fitness_ranking + pop[1].fitness
      table.insert(values_ranking, pop[1].fitness)
    end
    math.randomseed(s_Seeds[size])
    local gens_binary, fitness_binary, values_binary = 0, 0.0, {}
    for i = 1, RUNS do
      local gens, pop = EV(size, size, BinaryTournamentSelection)
      gens_binary = gens_binary + gens
      fitness_binary = fitness_binary + pop[1].fitness
      table.insert(values_binary, pop[1].fitness)
    end
    print(string.format("Pop Size: %d, Gens Infinite: %.2f, Gens Ranking/Binary Tournament: %.2f/%.2f, Fitness Ranking/Binary Tournament: %.2f/%.2f", size, gens_infinite, gens_ranking / RUNS, gens_binary / RUNS, fitness_ranking / RUNS, fitness_binary / RUNS))
    table.insert(points_finite, {x = size, y = gens_ranking / RUNS})
    table.insert(points_ranking, {x = size, y = gens_ranking})
    table.insert(points_binary, {x = size, y = gens_binary})
    table.insert(points_binary_fitness, {x = size, y = fitness_binary, values = ((#points_binary_fitness + 1) - points_binary_fitness.stddev_start) % points_binary_fitness.stddev_interval == 0 and values_binary})
    table.insert(points_ranking_fitness, {x = size, y = fitness_ranking, values = ((#points_ranking_fitness + 1) - points_ranking_fitness.stddev_start) % points_ranking_fitness.stddev_interval == 0 and values_ranking})
  end
  
  NormalizeGraphs(graphs3, RUNS)
  NormalizeGraphs(graphs4, RUNS)

  local bmp = Bitmap.new(IMAGE_WIDTH, IMAGE_HEIGHT, RGB_BLACK)
  DrawGraphs(bmp, graphs2, {int_x = true, skip_KP = true, min_x = 0, min_y = 0, div_x = 5, max_x = 250})
  bmp:WriteBMP(filename2)

  local bmp = Bitmap.new(IMAGE_WIDTH, IMAGE_HEIGHT, RGB_BLACK)
  DrawGraphs(bmp, graphs3, {int_x = true, skip_KP = true, min_x = 0, min_y = 0, div_x = 5, max_x = 250})
  bmp:WriteBMP(filename3)

  local bmp = Bitmap.new(IMAGE_WIDTH, IMAGE_HEIGHT, RGB_BLACK)
  DrawGraphs(bmp, graphs4, {int_x = true, skip_KP = true, min_x = 0, div_x = 5, div_y = 8, max_x = 250})
  bmp:WriteBMP(filename4)
end

local function TakeOverCurvesTruncationVSTournament(m, filename)
  local graphs = {name_x = string.format("Generation for m=%d", m), name_y = "TakeOver Proportion", funcs = {}}
  local name_trunc1 = string.format("Truncation k=1")
  local name_trunc5 = string.format("Truncation k=5")
  local name_tourn2 = string.format("Tournament q=2")
  local name_tourn3 = string.format("Tournament q=3")
  local name_tourn5 = string.format("Tournament q=5")
  local points_trunc1 = {color = RGB_CYAN}
  local points_trunc5 = {color = RGB_WHITE}
  local points_tourn2 = {color = RGB_RED}
  local points_tourn3 = {color = RGB_GREEN}
  local points_tourn5 = {color = RGB_BLUE}
  graphs.funcs[name_trunc1] = points_trunc1
  graphs.funcs[name_trunc5] = points_trunc5
  graphs.funcs[name_tourn2] = points_tourn2
  graphs.funcs[name_tourn3] = points_tourn3
  graphs.funcs[name_tourn5] = points_tourn5
  
  local Pt_trunc1, Pt_trunc5, Pt_tourn2, Pt_tourn3, Pt_tourn5 = 1.0 / m, 1.0 / m, 1.0 / m, 1.0 / m, 1.0 / m
  table.insert(points_trunc1, {x = 0, y = Pt_trunc1})
  table.insert(points_trunc5, {x = 0, y = Pt_trunc5})
  table.insert(points_tourn2, {x = 0, y = Pt_tourn2})
  table.insert(points_tourn3, {x = 0, y = Pt_tourn3})
  table.insert(points_tourn5, {x = 0, y = Pt_tourn5})
  for gen = 0, 12 do
    Pt_trunc1 = clamp((1.0 / m) * math.pow(m / 1, gen), 0, 1)
    Pt_trunc5 = clamp((1.0 / m) * math.pow(m / 5, gen), 0, 1)
    Pt_tourn2 = clamp(1.0 - math.pow(1.0 - Pt_tourn2, 2), 0, 1)
    Pt_tourn3 = clamp(1.0 - math.pow(1.0 - Pt_tourn3, 3), 0, 1)
    Pt_tourn5 = clamp(1.0 - math.pow(1.0 - Pt_tourn5, 5), 0, 1)
    table.insert(points_trunc1, {x = gen, y = Pt_trunc1})
    table.insert(points_trunc5, {x = gen, y = Pt_trunc5})
    table.insert(points_tourn2, {x = gen, y = Pt_tourn2})
    table.insert(points_tourn3, {x = gen, y = Pt_tourn3})
    table.insert(points_tourn5, {x = gen, y = Pt_tourn5})
  end
  
  local bmp = Bitmap.new(IMAGE_WIDTH, IMAGE_HEIGHT, RGB_BLACK)
  DrawGraphs(bmp, graphs, {int_x = true, skip_KP = true, min_x = 0, min_y = 0, div_x = 6, div_y = 5, max_y = 1.0})
  bmp:WriteBMP(filename)
end

local function CompareUniformTournamentFitnessProportional(filename, filename2)
  local graphs = {name_x = "Parent=Offspring Population Size", name_y = "Average Generations to Fixed-Point Convergence", funcs = {}}
  local name_uniform = string.format("Uniform Stochastic")
  local name_tournament = string.format("Binary Tournament")
  local name_fitness = string.format("Fitness-Proportional")
  local points_uniform = {color = RGB_GREEN}
  local points_tournament = {color = RGB_CYAN}
  local points_fitness = {color = RGB_MAGENTA}
  graphs.funcs[name_uniform] = points_uniform
  graphs.funcs[name_tournament] = points_tournament
  graphs.funcs[name_fitness] = points_fitness
  
  local graphs2 = {name_x = "Parent=Offspring Population Size", name_y = "Average Population Fitness at Fixed-Point Convergence", funcs = {}}
  local points2_uniform = {color = RGB_GREEN}
  local points2_tournament = {color = RGB_CYAN}
  local points2_fitness = {color = RGB_MAGENTA}
  graphs2.funcs[name_uniform] = points2_uniform
  graphs2.funcs[name_tournament] = points2_tournament
  graphs2.funcs[name_fitness] = points2_fitness
  
  for size = 10, MAX_POP_SIZE do
    math.randomseed(s_Seeds[size])
    local gens_uniform, fitness_uniform = 0, 0
    for run = 1, RUNS do
      local gens, pop = EV(size, size, UniformStochasticSelection)
      gens_uniform = gens_uniform + gens
      fitness_uniform = fitness_uniform + pop[1].fitness
    end
    table.insert(points_uniform, {x = size, y = gens_uniform})
    table.insert(points2_uniform, {x = size, y = fitness_uniform})
    
    math.randomseed(s_Seeds[size])
    local gens_tournament, fitness_tournament = 0, 0
    for run = 1, RUNS do
      local gens, pop = EV(size, size, BinaryTournamentSelection)
      gens_tournament = gens_tournament + gens
      fitness_tournament = fitness_tournament + pop[1].fitness
    end
    table.insert(points_tournament, {x = size, y = gens_tournament})
    table.insert(points2_tournament, {x = size, y = fitness_tournament})
    
    math.randomseed(s_Seeds[size])
    local gens_fitness, fitness_fitness = 0, 0
    for run = 1, RUNS do
      local gens, pop = EV(size, size, FitnessProportionalSelection)
      gens_fitness = gens_fitness + gens
      fitness_fitness = fitness_fitness + pop[1].fitness
    end
    table.insert(points_fitness, {x = size, y = gens_fitness})
    table.insert(points2_fitness, {x = size, y = fitness_fitness})

    print(string.format("Pop Size: %d, Gens Uniform: %.2f, Gens Tournament: %.2f, Gens Fitness-Proportional: %.2f", size, gens_uniform / RUNS, gens_tournament / RUNS, gens_fitness / RUNS))
  end
  
  NormalizeGraphs(graphs, RUNS)
  NormalizeGraphs(graphs2, RUNS)
  
  local bmp = Bitmap.new(IMAGE_WIDTH, IMAGE_HEIGHT, RGB_BLACK)
  DrawGraphs(bmp, graphs, {int_x = true, skip_KP = true, min_x = 0, min_y = 0, div_x = 5, div_y = 5})
  bmp:WriteBMP(filename)
  local bmp = Bitmap.new(IMAGE_WIDTH, IMAGE_HEIGHT, RGB_BLACK)
  DrawGraphs(bmp, graphs2, {int_x = true, skip_KP = true, min_x = 0, min_y = 0, div_x = 5, div_y = 5, max_y = MAX_FITNESS, y_line = MAX_FITNESS})
  bmp:WriteBMP(filename2)
end

local function CompareDeterministicVSStochastic(m, filename)
  local graphs = {name_x = "Offspring Population Size n", name_y = "Average Generations to Fixed-Point Convergence", funcs = {}}
  local name_deterministic = string.format("deterministic parent selection")
  local name_stochastic = string.format("stochastic parent selection")
  local points_stochastic = {color = RGB_GREEN}
  local points_deterministic = {color = RGB_CYAN}
  graphs.funcs[name_deterministic] = points_deterministic
  graphs.funcs[name_stochastic] = points_stochastic
  
  for brood_size = 1, 4 do
    math.randomseed(s_Seeds[brood_size])
    local gens_deterministic = 0
    if brood_size > 1 then
      for run = 1, RUNS do
        local gens = EV(m, brood_size * m, DeterministicSelection)
        gens_deterministic = gens_deterministic + gens
      end
    end
    table.insert(points_deterministic, {x = brood_size * m, y = gens_deterministic})
    
    math.randomseed(s_Seeds[brood_size])
    local gens_stochastic = 0
    for run = 1, RUNS do
      local gens = EV(m, brood_size * m, UniformStochasticSelection)
      gens_stochastic = gens_stochastic + gens
    end
    table.insert(points_stochastic, {x = brood_size * m, y = gens_stochastic})
    
    print(string.format("Offsping Size: %d, Gens Deterministic: %.2f, Gens Stochastic: %.2f", brood_size * m, gens_deterministic / RUNS, gens_stochastic / RUNS))
  end
  
  NormalizeGraphs(graphs, RUNS)
  
  local bmp = Bitmap.new(IMAGE_WIDTH, IMAGE_HEIGHT, RGB_BLACK)
  DrawGraphs(bmp, graphs, {int_x = true, skip_KP = true, min_x = 0, min_y = 0, div_x = 5, div_y = 4, max_x = 250, max_y = 200})
  bmp:WriteBMP(filename)
end

local function ParentUniformOffspringBinary(filename)
  local graphs = {name_x = "Offspring Population Size n", name_y = "Average Generations to Fixed-Point Convergence", funcs = {}}
  local m_values = {10, 20, 50}
  local colors = {RGB_GREEN, RGB_CYAN, RGB_RED}
  for i, m in ipairs(m_values) do
    local name = string.format("B-Tourn(m=%d,n=b*%d)", m, m);
    local points = {color = colors[i]}
    graphs.funcs[name] = points
    
    for brood_size = 1, MAX_POP_SIZE // m do
      local n = m * brood_size
      math.randomseed(s_Seeds[brood_size])
      local gens_total = 0
      for run = 1, RUNS do
        local gens = EV(m, n, UniformStochasticSelection, BinaryTournamentSelection)
        gens_total = gens_total + gens
      end
      table.insert(points, {x = n, y = gens_total})
      print(string.format("m=%d,n=%d, Gens: %.2f", m, n, gens_total / RUNS))
    end
  end
  
  local name = "B-Tourn(m=n)"
  local points = {color = RGB_WHITE}
  graphs.funcs[name] = points
  for m = 10, MAX_POP_SIZE do
    math.randomseed(s_Seeds[m])
    local gens_total = 0
    for run = 1, RUNS do
      local gens = EV(m, m, UniformStochasticSelection, BinaryTournamentSelection)
      gens_total = gens_total + gens
    end
    table.insert(points, {x = m, y = gens_total})
    print(string.format("m=%d, Gens: %.2f", m, gens_total / RUNS))
  end
  
  NormalizeGraphs(graphs, RUNS)
  
  local bmp = Bitmap.new(IMAGE_WIDTH, IMAGE_HEIGHT, RGB_BLACK)
  DrawGraphs(bmp, graphs, {int_x = true, skip_KP = true, min_x = 0, min_y = 0, div_x = 5, div_y = 4, max_x = 250, max_y = 20})
  bmp:WriteBMP(filename)
end

local function CompareOverlappingNaturalSelection(m, filename)
  local graphs = {name_x = "Offspring Population Size n", name_y = "Average Generations to Fixed-Point Convergence", funcs = {}}
  local name_deterministic = string.format("gen(%d,n) deterministic", m)
  local name_stochastic = string.format("gen(%d,n) stochaststic", m)
  local points_deterministic = {color = RGB_GREEN}
  local points_stochastic = {color = RGB_CYAN}
  graphs.funcs[name_deterministic] = points_deterministic
  graphs.funcs[name_stochastic] = points_stochastic
  
  for n = 1, MAX_POP_SIZE - m do
    local total_gens = 0
    for run = 1, RUNS do
      local gens = EV(m, n, UniformStochasticSelection, UniformStochasticSelection, "overlapping")
      total_gens = total_gens + gens
    end
    table.insert(points_stochastic, {x = n, y = total_gens})
    print(string.format("Stochastic+Stochastic: m=%d, n=%d, gens: %.2f", m, n, total_gens / RUNS))
  end
  
  for brood_size = 1, MAX_POP_SIZE // m do
    local n = m * brood_size
    local total_gens = 0
    for run = 1, RUNS do
      local gens = EV(m, n, DeterministicSelection, UniformStochasticSelection, "overlapping")
      total_gens = total_gens + gens
    end
    table.insert(points_deterministic, {x = n, y = total_gens})
    print(string.format("Deterministic+Stochastic: m=%d, n=%d, gens: %.2f", m, n, total_gens / RUNS))
  end
  
  NormalizeGraphs(graphs, RUNS)
  
  local bmp = Bitmap.new(IMAGE_WIDTH, IMAGE_HEIGHT, RGB_BLACK)
  DrawGraphs(bmp, graphs, {int_x = true, skip_KP = true, min_x = 0, min_y = 0, div_x = 5, div_y = 4, max_x = 250, max_y = 200})
  bmp:WriteBMP(filename)
end

local function CompareOverlappingBinaryTournament(m, filename)
  local graphs = {name_x = "Offspring Population Size", name_y = "Average Generations to Fixed-Point Convergence", funcs = {}}
  local name_btourn = string.format("B-Tourn(non-overlap) m=%d", m)
  local name_btourn_overlap = string.format("B-Tourn(overlap) m=%d", m)
  local name_btourn_mn = string.format("B-Tourn(non-overlap) m=n")
  local name_btourn_mn_overlap = string.format("B-Tourn(overlap) m=n")
  local points_btourn = {color = RGB_GREEN}
  local points_btourn_overlap = {color = RGB_CYAN}
  local points_btourn_mn = {color = RGB_RED}
  local points_btourn_mn_overlap = {color = RGB_WHITE}
  graphs.funcs[name_btourn] = points_btourn
  graphs.funcs[name_btourn_overlap] = points_btourn_overlap
  graphs.funcs[name_btourn_mn] = points_btourn_mn
  graphs.funcs[name_btourn_mn_overlap] = points_btourn_mn_overlap
  
  for n = 10, MAX_POP_SIZE do
    math.randomseed(s_Seeds[n])
    local total_gens = 0
    for run = 1, RUNS do
      local gens = EV(m, n, UniformStochasticSelection, BinaryTournamentSelection)
      total_gens = total_gens + gens
    end
    table.insert(points_btourn, {x = n, y = total_gens})
    print(string.format("Binary Tournament Non-Overlapping: m=%d, n=%d, gens: %.2f", m, n, total_gens / RUNS))
  end
  
  for n = 10, MAX_POP_SIZE do    
    math.randomseed(s_Seeds[n])
    local total_gens = 0
    for run = 1, RUNS do
      local gens = EV(m, n, UniformStochasticSelection, BinaryTournamentSelection, "overlapping")
      total_gens = total_gens + gens
    end
    table.insert(points_btourn_overlap, {x = n, y = total_gens})
    print(string.format("Binary Tournament Overlapping: m=%d, n=%d, gens: %.2f", m, n, total_gens / RUNS))
  end
  
  for m = 1, MAX_POP_SIZE do
    math.randomseed(s_Seeds[m])
    local total_gens = 0
    for run = 1, RUNS do
      local gens = EV(m, m, UniformStochasticSelection, BinaryTournamentSelection)
      total_gens = total_gens + gens
    end
    table.insert(points_btourn_mn, {x = m, y = total_gens})
    print(string.format("Binary Tournament Non-Overlapping: m=n=%d, gens: %.2f", m, total_gens / RUNS))
  end
  
  for m = 1, MAX_POP_SIZE do
    math.randomseed(s_Seeds[m])
    local total_gens = 0
    for run = 1, RUNS do
      local gens = EV(m, m, UniformStochasticSelection, BinaryTournamentSelection, "overlapping")
      total_gens = total_gens + gens
    end
    table.insert(points_btourn_mn_overlap, {x = m, y = total_gens})
    print(string.format("Binary Tournament Overlapping: m=n=%d, gens: %.2f", m, total_gens / RUNS))
  end
  
  NormalizeGraphs(graphs, RUNS)
  
  local bmp = Bitmap.new(IMAGE_WIDTH, IMAGE_HEIGHT, RGB_BLACK)
  DrawGraphs(bmp, graphs, {int_x = true, skip_KP = true, min_x = 0, min_y = 0, div_x = 5, div_y = 4, max_x = 250, max_y = 20})
  bmp:WriteBMP(filename)
end

local function CompareStandardEAsSelectionPressure(filename, filename2)
  local graphs = {name_x = "Offspring Population Size", name_y = "Average Generations to Fixed-Point Convergence", funcs = {}}
  local graphs2 = {name_x = "Offspring Population Size", name_y = "Average Population Fitness at Fixed-Point Convergence", funcs = {}}
  local name_ga_fprop = string.format("GA(Fitness-Proportional)")
  local name_ga_bt = string.format("GA(Binary Tournament)")
  local name_ep = string.format("EP")
  local name_es_plus = string.format("ES(plus)")
  local name_es_comma = string.format("ES(comma)")
  local points_ga_fprop, points2_ga_fprop = {color = RGB_GREEN}, {color = RGB_GREEN}
  local points_ga_bt, points2_ga_bt = {color = RGB_CYAN}, {color = RGB_CYAN}
  local points_ep, points2_ep = {color = RGB_WHITE}, {color = RGB_WHITE}
  local points_es_plus, points2_es_plus = {color = RGB_MAGENTA}, {color = RGB_MAGENTA}
  local points_es_comma, points2_es_comma = {color = RGB_ORANGE}, {color = RGB_ORANGE}
  graphs.funcs[name_ga_fprop], graphs2.funcs[name_ga_fprop] = points_ga_fprop, points2_ga_fprop
  graphs.funcs[name_ga_bt], graphs2.funcs[name_ga_bt] = points_ga_bt, points2_ga_bt
  graphs.funcs[name_ep], graphs2.funcs[name_ep] = points_ep, points2_ep
  graphs.funcs[name_es_plus], graphs2.funcs[name_es_plus] = points_es_plus, points2_es_plus
  graphs.funcs[name_es_comma], graphs2.funcs[name_es_comma] = points_es_comma, points2_es_comma
  
  for m = 10, MAX_POP_SIZE do
    math.randomseed(s_Seeds[m])
    local total_gens, total_fitness = 0, 0
    for run = 1, RUNS do
      local gens, pop = EV(m, 5 * m, UniformStochasticSelection, TruncateSelection)
      total_gens = total_gens + gens
      total_fitness = total_fitness + pop[1].fitness
    end
    table.insert(points_es_comma, {x = m, y = total_gens})
    table.insert(points2_es_comma, {x = m, y = total_fitness})
    print(string.format("ES(comma): m=%d, n=%d, gens: %.2f, fitness: %.2f", m, 5 * m, total_gens / RUNS, total_fitness / RUNS))
  end
  
  for m = 10, MAX_POP_SIZE do
    math.randomseed(s_Seeds[m])
    local total_gens, total_fitness = 0, 0
    for run = 1, RUNS do
      local gens, pop = EV(m, 5 * m, UniformStochasticSelection, TruncateSelection, "overlapping")
      total_gens = total_gens + gens
      total_fitness = total_fitness + pop[1].fitness
    end
    table.insert(points_es_plus, {x = m, y = total_gens})
    table.insert(points2_es_plus, {x = m, y = total_fitness})
    print(string.format("ES(plus): m=%d, n=%d, gens: %.2f, fitness: %.2f", m, 5 * m, total_gens / RUNS, total_fitness / RUNS))
  end
  
  for m = 10, MAX_POP_SIZE do
    math.randomseed(s_Seeds[m])
    local total_gens, total_fitness = 0, 0
    for run = 1, RUNS do
      local gens, pop = EV(m, m, DeterministicSelection, TruncateSelection, "overlapping")
      total_gens = total_gens + gens
      total_fitness = total_fitness + pop[1].fitness
    end
    table.insert(points_ep, {x = m, y = total_gens})
    table.insert(points2_ep, {x = m, y = total_fitness})
    print(string.format("EP: m=n=%d, gens: %.2f, fitness: %.2f", m, total_gens / RUNS, total_fitness / RUNS))
  end
  
  for m = 10, MAX_POP_SIZE do
    math.randomseed(s_Seeds[m])
    local total_gens, total_fitness = 0, 0
    for run = 1, RUNS do
      local gens, pop = GA(m, FitnessProportionalSelection)
      total_gens = total_gens + gens
      total_fitness = total_fitness + pop[1].fitness
    end
    table.insert(points_ga_fprop, {x = m, y = total_gens})
    table.insert(points2_ga_fprop, {x = m, y = total_fitness})
    print(string.format("GA Fitness-Proprtional: m=%d, gens: %.2f, fitness: %.2f", m, total_gens / RUNS, total_fitness / RUNS))
  end
  
  for m = 10, MAX_POP_SIZE do
    math.randomseed(s_Seeds[m])
    local total_gens, total_fitness = 0, 0
    for run = 1, RUNS do
      local gens, pop = GA(m, BinaryTournamentSelection)
      total_gens = total_gens + gens
      total_fitness = total_fitness + pop[1].fitness
    end
    table.insert(points_ga_bt, {x = m, y = total_gens})
    table.insert(points2_ga_bt, {x = m, y = total_fitness})
    print(string.format("GA Binary Tournament: m=%d, gens: %.2f, fitness: %.2f", m, total_gens / RUNS, total_fitness / RUNS))
  end
  
  NormalizeGraphs(graphs, RUNS)
  NormalizeGraphs(graphs2, RUNS)
  
  local bmp = Bitmap.new(IMAGE_WIDTH, IMAGE_HEIGHT, RGB_BLACK)
  DrawGraphs(bmp, graphs, {int_x = true, skip_KP = true, min_x = 0, min_y = 0, div_x = 5, div_y = 4, max_x = 250, max_y = 20})
  bmp:WriteBMP(filename)
  local bmp = Bitmap.new(IMAGE_WIDTH, IMAGE_HEIGHT, RGB_BLACK)
  DrawGraphs(bmp, graphs2, {int_x = true, skip_KP = true, min_x = 0, min_y = 92, div_x = 5, div_y = 8, max_x = 250, max_y = 20})
  bmp:WriteBMP(filename2)
end

local function CompareSUS(filename)
  local graphs = {name_x = "Parent Population Size", name_y = "Average Population Fitness at Fixed-Point Convergence", funcs = {}}
  local name_lr = string.format("Linear Ranking")
  local name_lr_sus = string.format("Linear Ranking SUS")
  local name_fp = string.format("Fitness-Proportional")
  local name_fp_sus = string.format("Fitness-Proportional SUS")
  local points_lr = {color = RGB_GREEN}
  local points_lr_sus = {color = RGB_MAGENTA}
  local points_fp = {color = RGB_CYAN}
  local points_fp_sus = {color = RGB_WHITE}
  graphs.funcs[name_lr] = points_lr
  graphs.funcs[name_lr_sus] = points_lr_sus
  graphs.funcs[name_fp] = points_fp
  graphs.funcs[name_fp_sus] = points_fp_sus
  
  for m = 10, MAX_POP_SIZE do
    math.randomseed(s_Seeds[m])
    local total_fitness = 0
    for run = 1, RUNS do
      local _, pop = EV(m, m, UniformStochasticSelection, FitnessProportionalSelection)
      total_fitness = total_fitness + pop[1].fitness
    end
    table.insert(points_fp, {x = m, y = total_fitness})
    print(string.format("Fitness-Proportional: m=%d, fitness: %.2f", m, total_fitness / RUNS))
  end
  
  for m = 10, MAX_POP_SIZE do
    math.randomseed(s_Seeds[m])
    local total_fitness = 0
    for run = 1, RUNS do
      local _, pop = EV(m, m, UniformStochasticSelection, SUS_FitnessProportionalSelection)
      total_fitness = total_fitness + pop[1].fitness
    end
    table.insert(points_fp_sus, {x = m, y = total_fitness})
    print(string.format("Fitness-Proportional SUS: m=%d, fitness: %.2f", m, total_fitness / RUNS))
  end
  
  for m = 10, MAX_POP_SIZE do
    math.randomseed(s_Seeds[m])
    local total_fitness = 0
    for run = 1, RUNS do
      local _, pop = EV(m, m, UniformStochasticSelection, LinearRankingSelection)
      total_fitness = total_fitness + pop[1].fitness
    end
    table.insert(points_lr, {x = m, y = total_fitness})
    print(string.format("Linear Ranking: m=%d, fitness: %.2f", m, total_fitness / RUNS))
  end
  
  for m = 10, MAX_POP_SIZE do
    math.randomseed(s_Seeds[m])
    local total_fitness = 0
    for run = 1, RUNS do
      local _, pop = EV(m, m, UniformStochasticSelection, SUS_LinearRanking)
      total_fitness = total_fitness + pop[1].fitness
    end
    table.insert(points_lr_sus, {x = m, y = total_fitness})
    print(string.format("Linear Ranking SUS: m=%d, fitness: %.2f", m, total_fitness / RUNS))
  end
  
  NormalizeGraphs(graphs, RUNS)
  
  local bmp = Bitmap.new(IMAGE_WIDTH, IMAGE_HEIGHT, RGB_BLACK)
  DrawGraphs(bmp, graphs, {int_x = true, skip_KP = true, min_x = 0, min_y = 90, div_x = 5, div_y = 4, max_x = 250})
  bmp:WriteBMP(filename)
end

FITNESS_FUNCTION = GetAssignedFitness
CompareNaturalSelection(IMAGE_FILENAME_FIGURE_6_1, IMAGE_FILENAME_FIGURE_6_2)
CompareTruncFactor(50, IMAGE_FILENAME_FIGURE_6_3, IMAGE_FILENAME_FIGURE_6_4)
TakeOverCurves(50, IMAGE_FILENAME_FIGURE_6_5, IMAGE_FILENAME_FIGURE_6_6, IMAGE_FILENAME_FIGURE_6_7, IMAGE_FILENAME_FIGURE_6_8)
TakeOverCurvesTruncationVSTournament(50, IMAGE_FILENAME_FIGURE_6_9)
CompareUniformTournamentFitnessProportional(IMAGE_FILENAME_FIGURE_6_10, IMAGE_FILENAME_FIGURE_6_11)
CompareDeterministicVSStochastic(50, IMAGE_FILENAME_FIGURE_6_12)
ParentUniformOffspringBinary(IMAGE_FILENAME_FIGURE_6_13)
CompareOverlappingNaturalSelection(50, IMAGE_FILENAME_FIGURE_6_14)
CompareOverlappingBinaryTournament(10, IMAGE_FILENAME_FIGURE_6_15)
CompareStandardEAsSelectionPressure(IMAGE_FILENAME_FIGURE_6_16, IMAGE_FILENAME_FIGURE_6_17)
CompareSUS(IMAGE_FILENAME_FIGURE_6_18)