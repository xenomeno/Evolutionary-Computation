dofile("Statistics.lua")
dofile("Bitmap.lua")
dofile("Graphics.lua")
dofile("EA_Common.lua")

local RUNS                          = 10
local BIRTHS                        = 1000
local BOUND_MIN                     = -5.0
local BOUND_MAX                     = 5.0
local BINARY_PRECISION              = 0.0001
local GENOME_BINARY_BITS            = math.ceil(math.log((BOUND_MAX - BOUND_MIN) / BINARY_PRECISION, 2))
local BINARY_MAX                    = math.pow(2, GENOME_BINARY_BITS) - 1
local CHROMOSOME_LENGTH             = 2 * GENOME_BINARY_BITS
local MUTATION_STEP                 = 0.1
local MUTATION_ADAPT_THRESHOLD      = 0.2     -- 20% or 1:5 rule
local PROB_CROSSOVER                = 0.6
local PROB_MUTATION                 = 0.005

local PRINT_DEBUG                   = false

local IMAGE_WIDTH                   = 1000
local IMAGE_HEIGHT                  = 1000
local IMAGE_FILENAME_FIGURE_3_1     = "EC/Figure3_1.bmp"
local IMAGE_FILENAME_FIGURE_3_2     = "EC/Figure3_2.bmp"
local IMAGE_FILENAME_FIGURE_3_3     = "EC/Figure3_3.bmp"
local IMAGE_FILENAME_FIGURE_3_4     = "EC/Figure3_4.bmp"
local IMAGE_FILENAME_FIGURE_3_5     = "EC/Figure3_5.bmp"
local IMAGE_FILENAME_FIGURE_3_6     = "EC/Figure3_6.bmp"
local IMAGE_FILENAME_FIGURE_3_7     = "EC/Figure3_7.bmp"
local IMAGE_FILENAME_FIGURE_3_8     = "EC/Figure3_8.bmp"
local IMAGE_FILENAME_FIGURE_3_9     = "EC/Figure3_9.bmp"
local IMAGE_FILENAME_FIGURE_3_10    = "EC/Figure3_10.bmp"

local function RealToBinary(x)
  return math.floor(BINARY_MAX * (x - BOUND_MIN) / (BOUND_MAX - BOUND_MIN))
end

local function BinaryToReal(x)
  return BOUND_MIN + (BOUND_MAX - BOUND_MIN) * x / BINARY_MAX
end

local function ChromToBitstring(chrom)
  local bits = {}
  for bit_pos = 1, CHROMOSOME_LENGTH do
    local bit_mask = 1 << (CHROMOSOME_LENGTH - bit_pos)
    bits[bit_pos] = ((chrom & bit_mask) ~= 0) and 1 or 0
  end
  
  return table.concat(bits, "")
end

local function Fitness(chrom, binary)
  local x1, x2
  if binary then
    x1, x2 = BinaryToReal((chrom >> GENOME_BINARY_BITS) & BINARY_MAX), BinaryToReal(chrom & BINARY_MAX)
  else
    x1, x2 = chrom.x1, chrom.x2
  end
  
  return x1 * x1 + x2 * x2
end

local function Copy(chrom)
  return {x1 = chrom.x1, x2 = chrom.x2}
end

local function Mutation(chrom)
  local x1 = clamp(chrom.x1 + GausMutation(MUTATION_STEP), BOUND_MIN, BOUND_MAX)
  local x2 = clamp(chrom.x2 + GausMutation(MUTATION_STEP), BOUND_MIN, BOUND_MAX)
  
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
end

local function GenInitPop(size, on_birth, binary)
  local pop = {binary = binary}
  while #pop < size do
    local chrom = {x1 = GenRand(BOUND_MIN, BOUND_MAX), x2 = GenRand(BOUND_MIN, BOUND_MAX)}
    if binary then
      chrom = (RealToBinary(chrom.x1) << GENOME_BINARY_BITS) + RealToBinary(chrom.x2)
    end
    table.insert(pop, {chrom = chrom, birthdate = #pop + 1})
    if on_birth then
      EvalPop(pop)
      on_birth(pop, #pop)
    end
  end

  return pop
end

local function EV(m, n, on_birth, mutation_step, mutation_adaptiveness)
  mutation_step = mutation_step or MUTATION_STEP
  
  local mutations_mul = 1.0
  local mutations_ratio = {[mutation_step] = {success = 0, total = 0}}
  
  local parents = GenInitPop(m, on_birth)
  local birth = #parents
  local gen = 1
  local offsprings = {}
  while birth < BIRTHS do
    local offsprings = {}
    for i = 1, n do
      local parent = parents[math.random(1, #parents)]
      local offspring = {chrom = Mutation(parent.chrom, mutation_step * mutations_mul)}
      offspring.fitness = Fitness(offspring.chrom)
      local ratio = mutations_ratio[mutation_step * mutations_mul]
      ratio.success = ratio.success + ((offspring.fitness > parent.fitness) and 1 or 0)
      ratio.total = ratio.total + 1
      table.insert(offsprings, offspring)
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
      if birth >= BIRTHS then break end
    end
    if mutation_adaptiveness then
      local ratio = mutations_ratio[mutation_step * mutations_mul]
      local success_mutations = ratio.success / ratio.total
      if success_mutations > MUTATION_ADAPT_THRESHOLD then
        mutations_mul = mutations_mul + mutation_adaptiveness
        mutations_ratio[mutation_step * mutations_mul] = mutations_ratio[mutation_step * mutations_mul] or {success = 0, total = 0}
      elseif success_mutations < MUTATION_ADAPT_THRESHOLD then
        mutations_mul = mutations_mul - mutation_adaptiveness
        mutations_ratio[mutation_step * mutations_mul] = mutations_ratio[mutation_step * mutations_mul] or {success = 0, total = 0}
      end
    end
    gen = gen + 1
  end
end

local function EP(m, on_birth, mutation_step)
  mutation_step = mutation_step or MUTATION_STEP
  
  local parents = GenInitPop(m, on_birth)
  table.sort(parents, function(ind1, ind2) return ind1.fitness > ind2.fitness end)
  local birth = #parents
  local gen = 1
  local offsprings = {}
  while birth < BIRTHS do
    local offsprings = {}
    for idx, parent in ipairs(parents) do
      local offspring = {chrom = Mutation(parent.chrom, mutation_step)}
      offspring.fitness = Fitness(offspring.chrom)
      offsprings[idx] = offspring
    end
    table.sort(offsprings, function(ind1, ind2) return ind1.fitness > ind2.fitness end)
    local new_parents = {}
    local parents_idx, offsprings_idx = 1, 1
    while #new_parents < #parents do
      if parents_idx > #parents then
        table.insert(new_parents, offsprings[offsprings_idx])
        offsprings_idx = offsprings_idx + 1
      elseif offsprings_idx > #offsprings then
        table.insert(new_parents, parents[parents_idx + 1])
      else
        local parent, offspring = parents[parents_idx], offsprings[offsprings_idx]
        if parent.fitness > offspring.fitness then
          table.insert(new_parents, parent)
          parents_idx = parents_idx + 1
        else
          table.insert(new_parents, offspring)
          offsprings_idx = offsprings_idx + 1
        end
      end
    end
    parents = new_parents
    EvalPop(parents)
    for k = 1, m do
      birth = birth + 1
      if on_birth then
        on_birth(parents, birth)
      end
    end
    gen = gen + 1
  end
end

local function ES(m, n, on_birth, mutation_step, mutation_adaptiveness)
  EV(m, n, on_birth, mutation_step, mutation_adaptiveness)
end

local function BinaryTournament(pop)
  local size = #pop
  local ind1 = pop[math.random(1, size)]
  local ind2 = pop[math.random(1, size)]
  
  return (ind1.fitness > ind2.fitness) and ind1 or ind2
end

local function Crossover(chrom1, chrom2)
  local offspring1_chrom = {x1 = chrom1.x1, x2 = chrom2.x2}
  local offspring2_chrom = {x1 = chrom2.x1, x2 = chrom1.x2}
  
  return {chrom = offspring1_chrom}, {chrom = offspring2_chrom}
end
  
local function GA(m, on_birth, mutation_step, prob_crossover)
  mutation_step = mutation_step or MUTATION_STEP
  
  local parents = GenInitPop(m, on_birth)
  local birth = #parents
  local gen = 1
  while birth < BIRTHS do
    local offsprings = {}
    while #offsprings < m do
      if prob_crossover then
        local parent1, parent2 = BinaryTournament(parents), BinaryTournament(parents)
        local offspring1, offspring2
        if math.random() < prob_crossover then
          offspring1, offspring2 = Crossover(parent1.chrom, parent2.chrom)
        else
          offspring1, offspring2 = {chrom = Copy(parent1.chrom)}, {chrom = Copy(parent2.chrom)}
        end
        offspring1.chrom, offspring2.chrom = Mutation(offspring1.chrom, mutation_step), Mutation(offspring2.chrom, mutation_step)
        table.insert(offsprings, offspring1)
        birth = birth + 1
        if on_birth then
          EvalPop(offsprings)
          on_birth(parents, birth)
        end
        if #offsprings < m then
          table.insert(offsprings, offspring2)
          birth = birth + 1
          if on_birth then
            EvalPop(offsprings)
            on_birth(parents, birth)
          end
        end
      else
        local parent = BinaryTournament(parents)
        local offspring = {chrom = Mutation(parent.chrom, mutation_step)}
        table.insert(offsprings, offspring)
        birth = birth + 1
        if on_birth then
          EvalPop(offsprings)
          on_birth(parents, birth)
        end
      end
    end
    parents = offsprings
    EvalPop(parents)
    gen = gen + 1
  end
end

local function BinaryCrossover(parent1, parent2)
  local xsite = math.random(1, CHROMOSOME_LENGTH - 1)
  local mask_full = (1 << CHROMOSOME_LENGTH) - 1
  local mask_right = (1 << (CHROMOSOME_LENGTH - xsite)) - 1
  local mask_left = mask_full - mask_right
  
  local offspring1 = (parent1 & mask_left) + (parent2 & mask_right)
  local offspring2 = (parent2 & mask_left) + (parent1 & mask_right)
  
  if PRINT_DEBUG then
    print(string.format("%s%d%s%d", string.rep(" ", xsite), xsite, string.rep(" ", CHROMOSOME_LENGTH - string.len(xsite) + 3), xsite))
    print(string.format("%s   %s", ChromToBitstring(parent1), ChromToBitstring(offspring1)))
    print(string.format("%s   %s", ChromToBitstring(parent2), ChromToBitstring(offspring2)))
  end
  
  return offspring1, offspring2
end


local function BinaryMutation(chrom, prob_mutation)
  for bit_pos = 1, CHROMOSOME_LENGTH do
    if math.random() < prob_mutation then
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

local function PrintPop(pop, gen)
  if not PRINT_DEBUG then return end
  
  print(string.format("Simulation time limit(# births): %d", BIRTHS))
  print(string.format("Using a genome with %d binary genes", GENOME_BINARY_BITS))
  print(string.format("Population size: %d", #pop))
  print("")
  print(string.format("Population data after %d births(generation %d)", #pop * gen, gen))
  print(string.format("Local fitness: Min: %.2f, Avg: %.2f, Max: %.2f", pop.min_fitness, pop.avg_fitness, pop.max_fitness))
  print(string.format("Individual  Birthdate  Fitness  %sGene Values", string.rep(" " , GENOME_BINARY_BITS)))
  for idx, ind in ipairs(pop) do
    local str = string.format("%8d    %7d    %3.2f", idx, ind.birthdate, ind.fitness)
    if pop.binary then
      print(string.format("%s%s%s", str, string.rep(" ", 32 - string.len(str)), ChromToBitstring(ind.chrom)))
    else
      print(string.format("%s%sx1=%.2f     x2=%.2f", str, string.rep(" ", 32 - string.len(str)), ind.chrom.x1, ind.chrom.x2))
    end
  end
  print("")
end

local function BinaryGA(m, on_birth, prob_crossover)
  local parents = GenInitPop(m, on_birth, "binary")
  local birth = #parents
  local gen = 1
  PrintPop(parents, gen)
  while birth < BIRTHS do
    local offsprings = {best_so_far = parents.best_so_far, binary = true}
    while #offsprings < m do
      if prob_crossover then
        local parent1, parent2 = BinaryTournament(parents), BinaryTournament(parents)
        local offspring1, offspring2
        if math.random() < prob_crossover then
          local chrom1, chrom2 = BinaryCrossover(parent1.chrom, parent2.chrom)
          offspring1, offspring2 = {chrom = chrom1}, {chrom = chrom2}
        else
          offspring1, offspring2 = {chrom = parent1.chrom}, {chrom = parent2.chrom}
        end
        offspring1.chrom, offspring2.chrom = BinaryMutation(offspring1.chrom, PROB_MUTATION), BinaryMutation(offspring2.chrom, PROB_MUTATION)
        table.insert(offsprings, offspring1)
        birth = birth + 1
        offspring1.birthdate = birth
        if on_birth then
          EvalPop(offsprings)
          on_birth(parents, birth)
        end
        if #offsprings < m then
          table.insert(offsprings, offspring2)
          birth = birth + 1
          offspring2.birthdate = birth
          if on_birth then
            EvalPop(offsprings)
            on_birth(parents, birth)
          end
        end
      else
        local parent = BinaryTournament(parents)
        local offspring = {chrom = BinaryMutation(parent.chrom, PROB_MUTATION), birthdate = birth + 1}
        table.insert(offsprings, offspring)
        birth = birth + 1
        if on_birth then
          EvalPop(offsprings)
          on_birth(parents, birth)
        end
      end
    end
    parents = offsprings
    EvalPop(parents)
    gen = gen + 1
    PrintPop(parents, gen)
  end
end

local function Figure3_1()
  local graphs =
  {
    funcs = {["EV(10, 1)"] = {color = RGB_CYAN}, ["EV(10, 10)"] = {color = RGB_GREEN}},
    name_x = "Number of Births",
    name_y = string.format("Fitness: Local Average(averaged over %d runs)", RUNS),
  }
  
  for run = 1, RUNS do
    EV(10, 1, function(pop, birth) PlotPop(pop, birth, graphs.funcs["EV(10, 1)"], "avg_fitness") end)
    EV(10, 10, function(pop, birth) PlotPop(pop, birth, graphs.funcs["EV(10, 10)"], "avg_fitness") end)
  end
  
  NormalizeGraphs(graphs, RUNS)
  
  local bmp = Bitmap.new(IMAGE_WIDTH, IMAGE_HEIGHT, RGB_BLACK)
  DrawGraphs(bmp, graphs, {int_x = true, skip_KP = true, min_y = 0.0})
  bmp:WriteBMP(IMAGE_FILENAME_FIGURE_3_1)
end

local function Figure3_2()
  local graphs =
  {
    funcs = {["EV(10, 10)"] = {color = RGB_CYAN}, ["EP(10)"] = {color = RGB_GREEN}},
    name_x = "Number of Births",
    name_y = string.format("Fitness: Local Average(averaged over %d runs)", RUNS),
  }
  
  for run = 1, RUNS do
    EV(10, 10, function(pop, birth) PlotPop(pop, birth, graphs.funcs["EV(10, 10)"], "avg_fitness") end)
    EP(10, function(pop, birth) PlotPop(pop, birth, graphs.funcs["EP(10)"], "avg_fitness") end)
  end
  
  NormalizeGraphs(graphs, RUNS)
  
  local bmp = Bitmap.new(IMAGE_WIDTH, IMAGE_HEIGHT, RGB_BLACK)
  DrawGraphs(bmp, graphs, {int_x = true, skip_KP = true, min_y = 0.0})
  bmp:WriteBMP(IMAGE_FILENAME_FIGURE_3_2)
end

local function Figure3_3()
  local graphs =
  {
    funcs = {["EV(10, 10)"] = {color = RGB_CYAN}, ["EP(10)"] = {color = RGB_GREEN}, ["ES(1, 10)"] = {color = RGB_RED}},
    name_x = "Number of Births",
    name_y = string.format("Fitness: Local Average(averaged over %d runs)", RUNS),
  }
  
  for run = 1, RUNS do
    ES(1, 10, function(pop, birth) PlotPop(pop, birth, graphs.funcs["ES(1, 10)"], "avg_fitness") end)
    EV(10, 10, function(pop, birth) PlotPop(pop, birth, graphs.funcs["EV(10, 10)"], "avg_fitness") end)
    EP(10, function(pop, birth) PlotPop(pop, birth, graphs.funcs["EP(10)"], "avg_fitness") end)
  end
  
  NormalizeGraphs(graphs, RUNS)
  
  local bmp = Bitmap.new(IMAGE_WIDTH, IMAGE_HEIGHT, RGB_BLACK)
  DrawGraphs(bmp, graphs, {int_x = true, skip_KP = true, min_y = 0.0})
  bmp:WriteBMP(IMAGE_FILENAME_FIGURE_3_3)
end

local function Figure3_4()
  local graphs =
  {
    funcs = {["EV(10, 10)"] = {color = RGB_CYAN}, ["EP(10)"] = {color = RGB_GREEN}, ["ES(5, 10)"] = {color = RGB_RED}},
    name_x = "Number of Births",
    name_y = string.format("Fitness: Local Average(averaged over %d runs)", RUNS),
  }
  
  for run = 1, RUNS do
    EV(10, 10, function(pop, birth) PlotPop(pop, birth, graphs.funcs["EV(10, 10)"], "avg_fitness") end)
    EP(10, function(pop, birth) PlotPop(pop, birth, graphs.funcs["EP(10)"], "avg_fitness") end)
    ES(5, 10, function(pop, birth) PlotPop(pop, birth, graphs.funcs["ES(5, 10)"], "avg_fitness") end)
  end
  
  NormalizeGraphs(graphs, RUNS)
  
  local bmp = Bitmap.new(IMAGE_WIDTH, IMAGE_HEIGHT, RGB_BLACK)
  DrawGraphs(bmp, graphs, {int_x = true, skip_KP = true, min_y = 0.0})
  bmp:WriteBMP(IMAGE_FILENAME_FIGURE_3_4)
end

local function Figure3_5()
  local graphs =
  {
    funcs = {["ES(1, 10, 0.1)"] = {color = RGB_RED}, ["ES(1, 10, 0.5)"] = {color = RGB_CYAN}, ["ES(1, 10, 5.0)"] = {color = RGB_GREEN}},
    name_x = "Number of Births",
    name_y = string.format("Fitness: Local Average(averaged over %d runs)", RUNS),
  }
  
  for run = 1, RUNS do
    ES(1, 10, function(pop, birth) PlotPop(pop, birth, graphs.funcs["ES(1, 10, 0.1)"], "avg_fitness") end, 0.1)
    ES(1, 10, function(pop, birth) PlotPop(pop, birth, graphs.funcs["ES(1, 10, 0.5)"], "avg_fitness") end, 0.5)
    ES(1, 10, function(pop, birth) PlotPop(pop, birth, graphs.funcs["ES(1, 10, 5.0)"], "avg_fitness") end, 5.0)
  end
  
  NormalizeGraphs(graphs, RUNS)
  
  local bmp = Bitmap.new(IMAGE_WIDTH, IMAGE_HEIGHT, RGB_BLACK)
  DrawGraphs(bmp, graphs, {int_x = true, skip_KP = true, min_y = 0.0})
  bmp:WriteBMP(IMAGE_FILENAME_FIGURE_3_5)
end

local function Figure3_6()
  local graphs =
  {
    funcs = {["ES(1, 10, 0.1)"] = {color = RGB_RED}, ["ES(1, 10, 0.5)"] = {color = RGB_CYAN}, ["ES(1, 10, 5.0)"] = {color = RGB_GREEN}, ["ES(1, 10, 1.0, 0.1)"] = {color = RGB_WHITE}},
    name_x = "Number of Births",
    name_y = string.format("Fitness: Local Average(averaged over %d runs)", RUNS),
  }
  
  for run = 1, RUNS do
    ES(1, 10, function(pop, birth) PlotPop(pop, birth, graphs.funcs["ES(1, 10, 0.1)"], "avg_fitness") end, 0.1)
    ES(1, 10, function(pop, birth) PlotPop(pop, birth, graphs.funcs["ES(1, 10, 0.5)"], "avg_fitness") end, 0.5)
    ES(1, 10, function(pop, birth) PlotPop(pop, birth, graphs.funcs["ES(1, 10, 5.0)"], "avg_fitness") end, 5.0)
    ES(1, 10, function(pop, birth) PlotPop(pop, birth, graphs.funcs["ES(1, 10, 1.0, 0.1)"], "avg_fitness") end, 1.0, 0.1)
  end
  
  NormalizeGraphs(graphs, RUNS)
  
  local bmp = Bitmap.new(IMAGE_WIDTH, IMAGE_HEIGHT, RGB_BLACK)
  DrawGraphs(bmp, graphs, {int_x = true, skip_KP = true, min_y = 0.0})
  bmp:WriteBMP(IMAGE_FILENAME_FIGURE_3_6)
end

local function Figure3_7()
  local graphs =
  {
    funcs = {["ES(1, 10)"] = {color = RGB_CYAN}, ["GA(10)"] = {color = RGB_GREEN}},
    name_x = "Number of Births",
    name_y = string.format("Fitness: Local Average(averaged over %d runs)", RUNS),
  }
  
  for run = 1, RUNS do
    ES(1, 10, function(pop, birth) PlotPop(pop, birth, graphs.funcs["ES(1, 10)"], "avg_fitness") end)
    GA(10, function(pop, birth) PlotPop(pop, birth, graphs.funcs["GA(10)"], "avg_fitness") end)
  end
  
  NormalizeGraphs(graphs, RUNS)
  
  local bmp = Bitmap.new(IMAGE_WIDTH, IMAGE_HEIGHT, RGB_BLACK)
  DrawGraphs(bmp, graphs, {int_x = true, skip_KP = true, min_y = 0.0})
  bmp:WriteBMP(IMAGE_FILENAME_FIGURE_3_7)
end

local function Figure3_8()
  local graphs =
  {
    funcs = {["ES(1, 10)"] = {color = RGB_CYAN}, ["GA(10)"] = {color = RGB_GREEN}, ["GA-X(10)"] = {color = RGB_WHITE}},
    name_x = "Number of Births",
    name_y = string.format("Fitness: Local Average(averaged over %d runs)", RUNS),
  }
  
  for run = 1, RUNS do
    ES(1, 10, function(pop, birth) PlotPop(pop, birth, graphs.funcs["ES(1, 10)"], "avg_fitness") end)
    GA(10, function(pop, birth) PlotPop(pop, birth, graphs.funcs["GA(10)"], "avg_fitness") end)
    GA(10, function(pop, birth) PlotPop(pop, birth, graphs.funcs["GA-X(10)"], "avg_fitness") end, nil, PROB_CROSSOVER)
  end
  
  NormalizeGraphs(graphs, RUNS)
  
  local bmp = Bitmap.new(IMAGE_WIDTH, IMAGE_HEIGHT, RGB_BLACK)
  DrawGraphs(bmp, graphs, {int_x = true, skip_KP = true, min_y = 0.0})
  bmp:WriteBMP(IMAGE_FILENAME_FIGURE_3_8)
end

local function Figure3_9()
  local graphs =
  {
    funcs = {["GA-B(10)"] = {color = RGB_GREEN}, ["GA-R(10)"] = {color = RGB_CYAN}},
    name_x = "Number of Births",
    name_y = string.format("Fitness: Local Average(averaged over %d runs)", RUNS),
  }
  
  for run = 1, RUNS do
    BinaryGA(10, function(pop, birth) PlotPop(pop, birth, graphs.funcs["GA-B(10)"], "avg_fitness") end, PROB_CROSSOVER)
    GA(10, function(pop, birth) PlotPop(pop, birth, graphs.funcs["GA-R(10)"], "avg_fitness") end, nil, PROB_CROSSOVER)
  end
  
  NormalizeGraphs(graphs, RUNS)
  
  local bmp = Bitmap.new(IMAGE_WIDTH, IMAGE_HEIGHT, RGB_BLACK)
  DrawGraphs(bmp, graphs, {int_x = true, skip_KP = true, min_y = 0.0})
  bmp:WriteBMP(IMAGE_FILENAME_FIGURE_3_9)
end

local function Figure3_10()
  local graphs =
  {
    funcs = {["GA-B(10)"] = {color = RGB_GREEN}, ["GA-R(10)"] = {color = RGB_CYAN}},
    name_x = "Number of Births",
    name_y = string.format("Fitness: Best so far(averaged over %d runs)", RUNS),
  }
  
  for run = 1, RUNS do 
    local max_bin_ga, max_ga
    BinaryGA(10, function(pop, birth)
      local points = graphs.funcs["GA-B(10)"]
      max_bin_ga = (not max_bin_ga or pop.max_fitness > max_bin_ga) and pop.max_fitness or max_bin_ga
      if points[birth] then
        points[birth].y = points[birth].y + max_bin_ga
      else
        points[birth] = {x = birth, y = max_bin_ga}
      end
    end, PROB_CROSSOVER)
    GA(10, function(pop, birth)
      local points = graphs.funcs["GA-R(10)"]
      max_ga = (not max_ga or pop.max_fitness > max_ga) and pop.max_fitness or max_ga
      if points[birth] then
        points[birth].y = points[birth].y + max_ga
      else
        points[birth] = {x = birth, y = max_ga}
      end
    end, nil, PROB_CROSSOVER)
  end
  
  NormalizeGraphs(graphs, RUNS)
  
  local bmp = Bitmap.new(IMAGE_WIDTH, IMAGE_HEIGHT, RGB_BLACK)
  DrawGraphs(bmp, graphs, {int_x = true, skip_KP = true, min_y = 0.0})
  bmp:WriteBMP(IMAGE_FILENAME_FIGURE_3_10)
end

Figure3_1()
Figure3_2()
Figure3_3()
Figure3_4()
Figure3_5()
Figure3_6()
Figure3_7()
Figure3_8()
Figure3_9()
Figure3_10()
