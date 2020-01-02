dofile("Statistics.lua")
dofile("Bitmap.lua")
dofile("Graphics.lua")
dofile("EC_Common.lua")

local RUNS                          = 10
local BIRTHS                        = 1000
local BOUND_MIN                     = -5.5
local BOUND_MAX                     = 4.5
local BINARY_PRECISION              = 0.0001
local GENOME_BINARY_BITS            = math.ceil(math.log((BOUND_MAX - BOUND_MIN) / BINARY_PRECISION, 2))
local BINARY_MAX                    = math.pow(2, GENOME_BINARY_BITS) - 1
local CHROMOSOME_LENGTH             = 2 * GENOME_BINARY_BITS
local MUTATION_STEP                 = 1.0
local PROB_CROSSOVER                = 1.0
local PROB_MUTATION                 = 1 / CHROMOSOME_LENGTH

local IMAGE_WIDTH                   = 1000
local IMAGE_HEIGHT                  = 1000
local IMAGE_FILENAME_FIGURE_4_2     = "EC/Figure4_2.bmp"
local IMAGE_FILENAME_FIGURE_4_3     = "EC/Figure4_3.bmp"
local IMAGE_FILENAME_FIGURE_4_4     = "EC/Figure4_4.bmp"
local IMAGE_FILENAME_FIGURE_4_6     = "EC/Figure4_6.bmp"
local IMAGE_FILENAME_FIGURE_4_7     = "EC/Figure4_7.bmp"
local IMAGE_FILENAME_FIGURE_4_8     = "EC/Figure4_8.bmp"
local IMAGE_FILENAME_FIGURE_4_10    = "EC/Figure4_10.bmp"
local IMAGE_FILENAME_FIGURE_4_11    = "EC/Figure4_11.bmp"

local function RealToBinary(x)
  return math.floor(BINARY_MAX * (x - BOUND_MIN) / (BOUND_MAX - BOUND_MIN))
end

local function BinaryToReal(x)
  return BOUND_MIN + (BOUND_MAX - BOUND_MIN) * x / BINARY_MAX
end

local function BinaryEncode(chrom)
  return (RealToBinary(chrom.x1) << GENOME_BINARY_BITS) + RealToBinary(chrom.x2)
end

local function BinaryDecode(chrom)
  return BinaryToReal((chrom >> GENOME_BINARY_BITS) & BINARY_MAX), BinaryToReal(chrom & BINARY_MAX)
end

local function Fitness(chrom, binary)
  local x1, x2
  if binary then
    x1, x2 = BinaryDecode(chrom)
  else
    x1, x2 = chrom.x1, chrom.x2
  end
  
  return x1 * x1 + x2 * x2
end

local function Copy(chrom)
  return {x1 = chrom.x1, x2 = chrom.x2}
end

local function Mutation(chrom, mutate_all)
  local x1, x2 = chrom.x2, chrom.x2
  if mutate_all or math.random() < 0.5 then
    x1 = clamp(x1 + GausMutation(MUTATION_STEP), BOUND_MIN, BOUND_MAX)
  end
  if mutate_all or math.random() < 0.5 then
    x2 = clamp(x2 + GausMutation(MUTATION_STEP), BOUND_MIN, BOUND_MAX)
  end
  
  return {x1 = x1, x2 = x2}
end

local function EvalPop(pop)
  local min_fitness, max_fitness
  local total_fitness = 0.0
  for _, ind in ipairs(pop) do
    local fitness = Fitness(ind.chrom, pop.binary)
    ind.fitness = fitness
    min_fitness = (not min_fitness or fitness < min_fitness) and fitness or min_fitness
    max_fitness = (not max_fitness or fitness > max_fitness) and fitness or max_fitness
    total_fitness = total_fitness + fitness
  end
  pop.min_fitness, pop.max_fitness, pop.avg_fitness = min_fitness, max_fitness, total_fitness / #pop
  pop.spread_fitness = max_fitness - min_fitness
end

local function GenInitPop(size, on_birth, binary)
  local pop = {binary = binary}
  while #pop < size do
    local chrom = {x1 = GenRand(BOUND_MIN, BOUND_MAX), x2 = GenRand(BOUND_MIN, BOUND_MAX)}
    if binary then
      chrom = BinaryEncode(chrom)
    end
    table.insert(pop, {chrom = chrom, birthdate = #pop + 1})
    if on_birth then
      EvalPop(pop)
      on_birth(pop, #pop)
    end
  end

  return pop
end

local function PickIndividual(prob, total_prob)
  local slot = math.random() * total_prob
  
  local left, right = 1, #prob
  while left + 1 < right do
    local middle = (left + right) // 2
    local total_prob = prob[middle].total_prob
    if slot == total_prob then
      return prob[middle].ind
    elseif slot < total_prob then
      right = middle
    else
      left = middle
    end
  end
  
  local idx = (slot < prob[left].total_prob + prob[left].ind_prob) and left or right
  
  return prob[idx].ind
end

local function PickPopulation(pool, size, total_prob)
  local pop = {}
  while #pop < size do
    table.insert(pop, PickIndividual(pool, total_prob))
  end
  
  return pop
end

local function UniformSelection(parents, offsprings, size)
  local prob = 1.0 / (#parents + #offsprings)
  
  local pool = {}
  local total_prob = 0.0
  for _, ind in ipairs(parents) do
    total_prob = total_prob + prob
    table.insert(pool, {ind_prob = prob, ind = ind, total_prob = total_prob})
  end
  for _, ind in ipairs(offsprings) do
    total_prob = total_prob + prob
    table.insert(pool, {ind_prob = prob, ind = ind, total_prob = total_prob})
  end
  
  return PickPopulation(pool, size, total_prob)
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

local function LinearRankingSelection(parents, offsprings, size)
  local pool = GetMergedSorted(parents, offsprings)
  
  local highest_prob = 0.1
  local total_prob = 0.0
  for idx, ind in ipairs(pool) do
    local ind_prob = highest_prob - highest_prob * (idx - 1) / (#pool - 1)
    total_prob = total_prob + ind_prob
    pool[idx] = {ind = ind, ind_prob = ind_prob, total_prob = total_prob}
  end
  
  return PickPopulation(pool, size, total_prob)
end

local function TruncateSelection(parents, offsprings, size)
  local pool = GetMergedSorted(parents, offsprings)
  local top = #pool // 5
  local prob = 1.0 / top
  local total_prob = 0.0
  for idx, ind in ipairs(pool) do
    local ind_prob = (idx <= top) and prob or 0.0
    total_prob = total_prob + ind_prob
    pool[idx] = {ind = ind, ind_prob = prob, total_prob = total_prob}
  end
  
  return PickPopulation(pool, size, total_prob)
end

local function EV(m, n, on_birth, selection, no_mutation, max_births, non_overlapping, mutale_all)
  max_births = max_births or BIRTHS
  
  local parents = GenInitPop(m, on_birth)
  local birth = #parents
  local gen = 1
  local offsprings = {}
  while birth < max_births do
    local offsprings = {}
    for i = 1, n do
      local parent = parents[math.random(1, #parents)]
      local chrom = no_mutation and Copy(parent.chrom) or Mutation(parent.chrom, mutate_all)
      local offspring = {chrom = chrom}
      offspring.fitness = Fitness(offspring.chrom)
      table.insert(offsprings, offspring)
    end
    if selection then
      parents = selection(non_overlapping and {} or parents, offsprings, #parents)
      EvalPop(parents)
      for i = 1, #parents do
        birth = birth + 1
        if on_birth then
          on_birth(parents, birth)
        end
        if birth >= max_births then break end
      end
    else
      for _, offspring in ipairs(offsprings) do
        local parent_idx = math.random(1, #parents)
        local parent = parents[parent_idx]
        parents[parent_idx] = (offspring.fitness > parent.fitness) and offspring or parent
        birth = birth + 1
        if on_birth then
          EvalPop(parents)
          on_birth(parents, birth)
        end
        if birth >= max_births then break end
      end
    end
    gen = gen + 1
  end
end

local function BinaryUniformCrossover(parent1, parent2, xcross)
  local offspring1, offspring2 = 0, 0
  local bit_mask = 1
  for bit_pos = 1, CHROMOSOME_LENGTH do
    local pos = CHROMOSOME_LENGTH - bit_pos + 1
    if xcross[pos] then
      offspring1 = offspring1 + (parent1 & bit_mask)
      offspring2 = offspring2 + (parent2 & bit_mask)
    else
      offspring1 = offspring1 + (parent2 & bit_mask)
      offspring2 = offspring2 + (parent1 & bit_mask)
    end
    bit_mask = bit_mask << 1
  end
  
  return offspring1, offspring2
end

local function BinaryParameterizedUniformCrossover(parent1, parent2, prob_head)
  prob_head = prob_head or 0.5
  
  local xcross = {}
  for pos = 1, CHROMOSOME_LENGTH do
    xcross[pos] = math.random() < prob_head
  end
  
  return BinaryUniformCrossover(parent1, parent2, xcross)
end

local function BinaryOnePointCrossover(parent1, parent2)
  local xsite = math.random(1, CHROMOSOME_LENGTH)
  local xcross = {}
  for i = 1, CHROMOSOME_LENGTH do
    xcross[i] = i < xsite
  end
  
  return BinaryUniformCrossover(parent1, parent2, xcross)
end

local function BinaryTwoPointsCrossover(parent1, parent2)
  local xsite1 = math.random(1, CHROMOSOME_LENGTH - 1)
  local xsite2 = math.random(xsite1 + 1, CHROMOSOME_LENGTH)
  local xcross = {}
  for i = 1, CHROMOSOME_LENGTH do
    xcross[i] = (i < xsite1) or (i >= xsite2)
  end
  
  return BinaryUniformCrossover(parent1, parent2, xcross)
end

local function BinaryMutation(chrom)
  for bit_pos = 1, CHROMOSOME_LENGTH do
    if math.random() < PROB_MUTATION then
      local bit_mask = 1 << (CHROMOSOME_LENGTH - bit_pos)
      if (chrom & bit_mask) ~= 0 then
        chrom = chrom - bit_mask
      else
        chrom = chrom + bit_mask
      end
    end
  end
  
  return chrom
end

local function BinaryEV(m, n, on_birth, crossover)
  max_births = max_births or BIRTHS
  
  local parents = GenInitPop(m, on_birth, "binary")
  local birth = #parents
  local gen = 1
  local offsprings = {}
  while birth < max_births do
    local offsprings = {}
    while #offsprings < n do
      local parent1 = parents[math.random(1, #parents)]
      local parent2 = parents[math.random(1, #parents)]
      local offspring1, offspring2
      if math.random() < PROB_CROSSOVER then
        offspring1, offspring2 = crossover(parent1.chrom, parent2.chrom)
      else
        offspring1, offspring2 = parent1.chrom, parent2.chrom
      end
      offspring1 = BinaryMutation(offspring1)
      table.insert(offsprings, {chrom = offspring1, fitness = Fitness(offspring1, "binary")})
      if #offsprings < n then
        offspring2 = BinaryMutation(offspring2)
        table.insert(offsprings, {chrom = offspring2, fitness = Fitness(offspring2, "binary")})
      end
    end
    for _, offspring in ipairs(offsprings) do
      local parent_idx = math.random(1, #parents)
      local parent = parents[parent_idx]
      parents[parent_idx] = (offspring.fitness > parent.fitness) and offspring or parent
      birth = birth + 1
      if on_birth then
        EvalPop(parents)
        on_birth(parents, birth)
      end
      if birth >= max_births then break end
    end
    gen = gen + 1
  end
end

local s_Seeds = {}
for run = 1, RUNS do
  s_Seeds[run] = math.random(1, 100000)
end
  
local function Figure4_2()
  local name10, name40, name100 = "EV(10, 1, 1.0)", "EV(40, 1, 1.0)", "EV(100, 1, 1.0)"
  local graphs =
  {
    funcs = {[name10] = {color = RGB_CYAN}, [name40] = {color = RGB_GREEN}, [name100] = {color = RGB_RED}},
    name_x = "Number of Births",
    name_y = string.format("Fitness: Best-So-Far(averaged over %d runs)", RUNS),
  }
  
  for run = 1, RUNS do
    local best_so_far = {}
    math.randomseed(s_Seeds[run])
    EV(10, 1, function(pop, birth)
      local max = pop.max_fitness
      best_so_far[name10] = (not best_so_far[name10] or max > best_so_far[name10]) and max or best_so_far[name10]
      points = graphs.funcs[name10]
      if points[birth] then
        points[birth].y = points[birth].y + max
      else
        points[birth] = {x = birth, y = max}
      end
    end)
    math.randomseed(s_Seeds[run])
    EV(40, 1, function(pop, birth)
      local max = pop.max_fitness
      best_so_far[name40] = (not best_so_far[name40] or max > best_so_far[name40]) and max or best_so_far[name40]
      points = graphs.funcs[name40]
      if points[birth] then
        points[birth].y = points[birth].y + max
      else
        points[birth] = {x = birth, y = max}
      end
    end)
    math.randomseed(s_Seeds[run])
    EV(100, 1, function(pop, birth)
      local max = pop.max_fitness
      best_so_far[name100] = (not best_so_far[name100] or max > best_so_far[name100]) and max or best_so_far[name100]
      points = graphs.funcs[name100]
      if points[birth] then
        points[birth].y = points[birth].y + max
      else
        points[birth] = {x = birth, y = max}
      end
    end)
  end
  
  NormalizeGraphs(graphs, RUNS)
  
  local bmp = Bitmap.new(IMAGE_WIDTH, IMAGE_HEIGHT, RGB_BLACK)
  DrawGraphs(bmp, graphs, {int_x = true, skip_KP = true, min_y = 0.0})
  bmp:WriteBMP(IMAGE_FILENAME_FIGURE_4_2)
end

local function Figure4_3()
  local name10, name40 = "EV(10, 1, 1.0)", "EV(40, 1, 1.0)"
  local graphs =
  {
    funcs = {[name10] = {color = RGB_CYAN, stddev_interval = 50, stddev_start = 20}, [name40] = {color = RGB_GREEN, stddev_interval = 50, stddev_start = 1}},
    name_x = "Number of Births",
    name_y = string.format("Fitness: Best-So-Far(averaged over %d runs)", RUNS),
  }
  
  for run = 1, RUNS do
    local best_so_far = {}
    math.randomseed(s_Seeds[run])
    EV(10, 1, function(pop, birth)
      local max = pop.max_fitness
      best_so_far[name10] = (not best_so_far[name10] or max > best_so_far[name10]) and max or best_so_far[name10]
      points = graphs.funcs[name10]
      if points[birth] then
        points[birth].y = points[birth].y + max
        table.insert(points[birth].values, max)
      else
        points[birth] = {x = birth, y = max, values = {max}}
      end
    end)
    math.randomseed(s_Seeds[run])
    EV(40, 1, function(pop, birth)
      local max = pop.max_fitness
      best_so_far[name40] = (not best_so_far[name40] or max > best_so_far[name40]) and max or best_so_far[name40]
      points = graphs.funcs[name40]
      if points[birth] then
        points[birth].y = points[birth].y + max
        table.insert(points[birth].values, max)
      else
        points[birth] = {x = birth, y = max, values = {max}}
      end
    end)
  end
  
  NormalizeGraphs(graphs, RUNS)
  
  local bmp = Bitmap.new(IMAGE_WIDTH, IMAGE_HEIGHT, RGB_BLACK)
  DrawGraphs(bmp, graphs, {int_x = true, skip_KP = true, min_y = 0.0})
  bmp:WriteBMP(IMAGE_FILENAME_FIGURE_4_3)
end

local function Figure4_4()
  local name10_1, name10_10, name10_20, name40_1 = "EV(10, 1, 1.0)", "EV(10, 10, 1.0)", "EV(10, 20, 1.0)", "EV(40, 1, 1.0)"
  local graphs =
  {
    funcs = {[name10_1] = {color = RGB_CYAN}, [name10_10] = {color = RGB_GREEN}, [name10_20] = {color = RGB_RED}, [name40_1] = {color = RGB_WHITE}},
    name_x = "Number of Births",
    name_y = string.format("Fitness: Best-So-Far(averaged over %d runs)", RUNS),
  }
  
  for run = 1, RUNS do
    local best_so_far = {}
    math.randomseed(s_Seeds[run])
    EV(10, 1, function(pop, birth)
      local max = pop.max_fitness
      best_so_far[name10_1] = (not best_so_far[name10_1] or max > best_so_far[name10_1]) and max or best_so_far[name10_1]
      points = graphs.funcs[name10_1]
      if points[birth] then
        points[birth].y = points[birth].y + max
      else
        points[birth] = {x = birth, y = max}
      end
    end)
    math.randomseed(s_Seeds[run])
    EV(10, 10, function(pop, birth)
      local max = pop.max_fitness
      best_so_far[name10_10] = (not best_so_far[name10_10] or max > best_so_far[name10_10]) and max or best_so_far[name10_10]
      points = graphs.funcs[name10_10]
      if points[birth] then
        points[birth].y = points[birth].y + max
      else
        points[birth] = {x = birth, y = max}
      end
    end)
    math.randomseed(s_Seeds[run])
    EV(10, 20, function(pop, birth)
      local max = pop.max_fitness
      best_so_far[name10_20] = (not best_so_far[name10_20] or max > best_so_far[name10_20]) and max or best_so_far[name10_20]
      points = graphs.funcs[name10_20]
      if points[birth] then
        points[birth].y = points[birth].y + max
      else
        points[birth] = {x = birth, y = max}
      end
    end)
    EV(40, 1, function(pop, birth)
      local max = pop.max_fitness
      best_so_far[name40_1] = (not best_so_far[name40_1] or max > best_so_far[name40_1]) and max or best_so_far[name40_1]
      points = graphs.funcs[name40_1]
      if points[birth] then
        points[birth].y = points[birth].y + max
      else
        points[birth] = {x = birth, y = max}
      end
    end)
  end
  
  NormalizeGraphs(graphs, RUNS)
  
  local bmp = Bitmap.new(IMAGE_WIDTH, IMAGE_HEIGHT, RGB_BLACK)
  DrawGraphs(bmp, graphs, {int_x = true, skip_KP = true, center_y = 0.0, max_y = 65.0})
  bmp:WriteBMP(IMAGE_FILENAME_FIGURE_4_4)
end

local function Figure4_6_And_4_7()
  local trunk, rank, uniform = "EV(10,10,Truncate Selection)", "EV(10,10,Linear Ranking Selection)", "EV(10,10,Uniform Selection)"
  local average =
  {
    funcs = {[trunk] = {color = RGB_CYAN}, [rank] = {color = RGB_GREEN}, [uniform] = {color = RGB_RED}},
    name_x = "Number of Births",
    name_y = string.format("Average Fitness(averaged over %d runs)", RUNS),
  }
  local spread =
  {
    funcs = {[trunk] = {color = RGB_CYAN}, [rank] = {color = RGB_GREEN}, [uniform] = {color = RGB_RED}},
    name_x = "Number of Births",
    name_y = string.format("Fitness Spread Max-Min(averaged over %d runs)", RUNS),
  }
  
  for run = 1, RUNS do
    math.randomseed(s_Seeds[run])
    EV(10, 10, function(pop, birth)
      PlotPop(pop, birth, average.funcs[uniform], "avg_fitness")
      PlotPop(pop, birth, spread.funcs[uniform], "spread_fitness")
    end, UniformSelection, "no mutation", 500)
    math.randomseed(s_Seeds[run])
    EV(10, 10, function(pop, birth)
      PlotPop(pop, birth, average.funcs[trunk], "avg_fitness")
      PlotPop(pop, birth, spread.funcs[trunk], "spread_fitness")
    end, TruncateSelection, "no mutation", 500)
    math.randomseed(s_Seeds[run])
    EV(10, 10, function(pop, birth)
      PlotPop(pop, birth, average.funcs[rank], "avg_fitness")
      PlotPop(pop, birth, spread.funcs[rank], "spread_fitness")
    end, LinearRankingSelection, "no mutation", 500)
  end
  
  NormalizeGraphs(average, RUNS)
  NormalizeGraphs(spread, RUNS)
  
  local bmp = Bitmap.new(IMAGE_WIDTH, IMAGE_HEIGHT, RGB_BLACK)
  DrawGraphs(bmp, average, {int_x = true, skip_KP = true, center_y = 0.0, max_y = 40.0})
  bmp:WriteBMP(IMAGE_FILENAME_FIGURE_4_6)
  local bmp = Bitmap.new(IMAGE_WIDTH, IMAGE_HEIGHT, RGB_BLACK)
  DrawGraphs(bmp, spread, {int_x = true, skip_KP = true, center_y = 0.0, max_y = 40.0})
  bmp:WriteBMP(IMAGE_FILENAME_FIGURE_4_7)
end

local function Figure4_8()
  local overlapping, non_overlapping = "EV(10,20,Overlapping)", "EV(10,20,Non-Overlapping)"
  local graphs =
  {
    funcs = {[overlapping] = {color = RGB_CYAN}, [non_overlapping] = {color = RGB_GREEN}},
    name_x = "Number of Births",
    name_y = string.format("Average Fitness(averaged over %d runs)", RUNS),
  }
  
  for run = 1, RUNS do
    math.randomseed(s_Seeds[run])
    EV(10, 20, function(pop, birth)
      PlotPop(pop, birth, graphs.funcs[overlapping], "avg_fitness")
    end, LinearRankingSelection, "no mutation", 300)
    math.randomseed(s_Seeds[run])
    EV(10, 20, function(pop, birth)
      PlotPop(pop, birth, graphs.funcs[non_overlapping], "avg_fitness")
    end, LinearRankingSelection, "no mutation", 300, "non-overlapping")
  end
  
  NormalizeGraphs(graphs, RUNS)
  
  local bmp = Bitmap.new(IMAGE_WIDTH, IMAGE_HEIGHT, RGB_BLACK)
  DrawGraphs(bmp, graphs, {int_x = true, skip_KP = true, center_y = 0.0, max_y = 40.0})
  bmp:WriteBMP(IMAGE_FILENAME_FIGURE_4_8)
end

local function Figure4_10()
  local mutate_all, mutate_one = "EV(10,10,Mutate All)", "EV(10,10,Mutate One)"
  local graphs =
  {
    funcs = {[mutate_all] = {color = RGB_CYAN}, [mutate_one] = {color = RGB_GREEN}},
    name_x = "Number of Births",
    name_y = string.format("Best-So-Far Fitness(averaged over %d runs)", RUNS),
  }
  
  for run = 1, RUNS do
    local best_mutate_all
    math.randomseed(s_Seeds[run])
    EV(10, 10, function(pop, birth)
      best_mutate_all = (not best_mutate_all or pop.max_fitness > best_mutate_all) and pop.max_fitness or best_mutate_all
      local points = graphs.funcs[mutate_all]
      if points[birth] then
        points[birth].y = points[birth].y + best_mutate_all
      else
        points[birth] = {x = birth, y = best_mutate_all}
      end
    end, nil, nil, nil, nil, "mutate all")
    local best_mutate_one
    math.randomseed(s_Seeds[run])
    EV(10, 10, function(pop, birth)
      best_mutate_one = (not best_mutate_one or pop.max_fitness > best_mutate_one) and pop.max_fitness or best_mutate_one
      local points = graphs.funcs[mutate_one]
      if points[birth] then
        points[birth].y = points[birth].y + best_mutate_one
      else
        points[birth] = {x = birth, y = best_mutate_one}
      end
    end)
  end
  
  NormalizeGraphs(graphs, RUNS)
  
  local bmp = Bitmap.new(IMAGE_WIDTH, IMAGE_HEIGHT, RGB_BLACK)
  DrawGraphs(bmp, graphs, {int_x = true, skip_KP = true, center_y = 0.0, max_y = 40.0})
  bmp:WriteBMP(IMAGE_FILENAME_FIGURE_4_10)
end

local function Figure4_11()
  local p02, p05, cp1, cp2 = "EV(10,10,0.2)", "EV(10,10,0.5)", "EV(10,10,One-Point Crossover)", "EV(10,10,Two-Point Crossover"
  local graphs =
  {
    funcs = {[p02] = {color = RGB_CYAN}, [p05] = {color = RGB_GREEN}, [cp1] = {color = RGB_RED}, [cp2] = {color = RGB_WHITE}},
    name_x = "Number of Births",
    name_y = string.format("Best-So-Far Fitness(averaged over %d runs)", RUNS),
  }
  
  for run = 1, RUNS do
    local best_so_far = {}
    math.randomseed(s_Seeds[run])
    BinaryEV(10, 10, function(pop, birth)
      best_so_far[p02] = (not best_so_far[p02] or pop.max_fitness > best_so_far[p02]) and pop.max_fitness or best_so_far[p02]
      local points = graphs.funcs[p02]
      if points[birth] then
        points[birth].y = points[birth].y + best_so_far[p02]
      else
        points[birth] = {x = birth, y = best_so_far[p02]}
      end
    end, function(parent1, parent2) return BinaryParameterizedUniformCrossover(parent1, parent2, 0.2) end)
    math.randomseed(s_Seeds[run])
    BinaryEV(10, 10, function(pop, birth)
      best_so_far[p05] = (not best_so_far[p05] or pop.max_fitness > best_so_far[p05]) and pop.max_fitness or best_so_far[p05]
      local points = graphs.funcs[p05]
      if points[birth] then
        points[birth].y = points[birth].y + best_so_far[p05]
      else
        points[birth] = {x = birth, y = best_so_far[p05]}
      end
    end, function(parent1, parent2) return BinaryParameterizedUniformCrossover(parent1, parent2, 0.5) end)
    math.randomseed(s_Seeds[run])
    BinaryEV(10, 10, function(pop, birth)
      best_so_far[cp1] = (not best_so_far[cp1] or pop.max_fitness > best_so_far[cp1]) and pop.max_fitness or best_so_far[cp1]
      local points = graphs.funcs[cp1]
      if points[birth] then
        points[birth].y = points[birth].y + best_so_far[cp1]
      else
        points[birth] = {x = birth, y = best_so_far[cp1]}
      end
    end, BinaryOnePointCrossover)
    math.randomseed(s_Seeds[run])
    BinaryEV(10, 10, function(pop, birth)
      best_so_far[cp2] = (not best_so_far[cp2] or pop.max_fitness > best_so_far[cp2]) and pop.max_fitness or best_so_far[cp2]
      local points = graphs.funcs[cp2]
      if points[birth] then
        points[birth].y = points[birth].y + best_so_far[cp2]
      else
        points[birth] = {x = birth, y = best_so_far[cp2]}
      end
    end, BinaryTwoPointsCrossover)
  end
  
  NormalizeGraphs(graphs, RUNS, 40.0)
  
  local bmp = Bitmap.new(IMAGE_WIDTH, IMAGE_HEIGHT, RGB_BLACK)
  DrawGraphs(bmp, graphs, {int_x = true, skip_KP = true, center_y = 0.0, center_y = 40.0})
  bmp:WriteBMP(IMAGE_FILENAME_FIGURE_4_11)
end

Figure4_2()
Figure4_3()
Figure4_4()
Figure4_6_And_4_7()
Figure4_8()
Figure4_10()
Figure4_11()