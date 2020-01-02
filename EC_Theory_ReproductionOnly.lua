dofile("Statistics.lua")
dofile("Bitmap.lua")
dofile("Graphics.lua")
dofile("EC_Common.lua")
dofile("GA_Common.lua")

local RUNS                          = 10
local POPULATION_SIZE               = 128
local MAX_GENERATIONS               = 100

local GENOME_BITS                   = 5
local GENOME_CARDINALITY            = math.pow(2, GENOME_BITS)
local GENOME_UNIFORM_FREQUENCY      = 1.0 / GENOME_CARDINALITY

local GENOME_REALS                  = 10
local GENOME_REAL_MIN               = -5.0
local GENOME_REAL_MAX               = 5.0

local IMAGE_WIDTH                   = 1000
local IMAGE_HEIGHT                  = 1000
local IMAGE_FILENAME_FIGURE_6_19    = "EC/Figure6_19.bmp"
local IMAGE_FILENAME_FIGURE_6_20    = "EC/Figure6_20.bmp"
local IMAGE_FILENAME_FIGURE_6_21    = "EC/Figure6_21.bmp"
local IMAGE_FILENAME_FIGURE_6_22    = "EC/Figure6_22.bmp"
local IMAGE_FILENAME_FIGURE_6_23    = "EC/Figure6_23.bmp"
local IMAGE_FILENAME_FIGURE_6_24    = "EC/Figure6_24.bmp"
local IMAGE_FILENAME_FIGURE_6_25    = "EC/Figure6_25.bmp"
local IMAGE_FILENAME_FIGURE_6_26    = "EC/Figure6_26.bmp"
local IMAGE_FILENAME_FIGURE_6_27    = "EC/Figure6_27.bmp"
local IMAGE_FILENAME_FIGURE_6_28    = "EC/Figure6_28.bmp"
local IMAGE_FILENAME_FIGURE_6_29    = "EC/Figure6_29.bmp"

local function RMS(pop)
  local c = {}
  for i = 0, GENOME_CARDINALITY - 1 do
    c[i] = 0
  end
  for _, ind in ipairs(pop) do
    local pow2 = GENOME_CARDINALITY >> 1
    local genotype = 0
    for _, allele in ipairs(ind) do
      genotype = genotype + ((allele == 'B') and pow2 or 0)
      pow2 = pow2 >> 1
    end
    c[genotype] = c[genotype] + 1
  end

  local RMS2 = 0
  for genotype = 0, GENOME_CARDINALITY - 1 do
    local frequency = c[genotype] / POPULATION_SIZE
    RMS2 = RMS2 + (frequency - GENOME_UNIFORM_FREQUENCY) * (frequency - GENOME_UNIFORM_FREQUENCY)
  end
  
  return math.sqrt(RMS2 / GENOME_CARDINALITY)
end

local function RMS_Hamming(pop)
  local rms, unique_pairs = 0, 0
  for i, genotype in ipairs(pop) do
    for j = i + 1, #pop do
      local genotype2 = pop[j]
      local dist = 0
      for k, allele in ipairs(genotype) do
        dist = dist + ((allele ~= genotype2[k]) and 1 or 0)
      end
      if dist > 0 then
        rms = rms + dist
        unique_pairs = unique_pairs + 1
      end
    end
  end
  
  return (unique_pairs == 0) and 0 or (rms / unique_pairs)
end

local function Dispersion(pop)
  local centroid = {}
  for _, genotype in ipairs(pop) do
    for i, allele in ipairs(genotype) do
      centroid[i] = (centroid[i] or 0) + allele
    end
  end
  local inv = 1.0 / #pop
  for i, allele in ipairs(centroid) do
    centroid[i] = allele * inv
  end
  
  local total_dist = 0
  for _, genotype in ipairs(pop) do
    local dist2 = 0
    for i, allele in ipairs(genotype) do
      dist2 = dist2 + (allele - centroid[i]) * (allele - centroid[i])
    end
    total_dist = total_dist + math.sqrt(dist2)
  end
  
  return total_dist * inv
end
  
local function GenInitEdgePop(size)
  local ind_A, ind_B = {}, {}
  for i = 1, GENOME_BITS do
    ind_A[i], ind_B[i] = 'A', 'B'
  end
  
  local pop = {}
  for i = 1, size // 2 do
    table.insert(pop, table.copy(ind_A))
    table.insert(pop, table.copy(ind_B))
  end

  return pop
end

local function GenInitUniformPop(size)
  local ind_A = {}
  for i = 1, GENOME_BITS do
    ind_A[i] = 'A'
  end
  
  local pop = {}
  for i = 1, size do
    table.insert(pop, table.copy(ind_A))
  end
  
  return pop
end

local function GenInitRandomPop(size)
  local pop = {}
  while #pop < size do
    local chrom = {}
    for i = 1, GENOME_BITS do
      chrom[i] = (math.random() < 0.5) and 'A' or 'B'
    end
    table.insert(pop, chrom)
  end
  
  return pop
end

local function GenInitRealRandomPop(size)
  local pop = {}
  while #pop < size do
    local genotype = {}
    for i = 1, GENOME_REALS do
      genotype[i] = GENOME_REAL_MIN + math.random() * (GENOME_REAL_MAX - GENOME_REAL_MIN)
    end
    table.insert(pop, genotype)
  end

  return pop
end

local function UniformStochasticSelection(parents, size)
  local offsprings = {}
  for i = 1, size do
    local parent = parents[math.random(1, #parents)]
    table.insert(offsprings, parent)
  end
  
  return offsprings
end

local function GenPermutation(size)
  local permutation = {}
  for i = 1, size do
    permutation[i] = i
  end
  for i = 1, #permutation do
    local j = math.random(1, #permutation)
    permutation[i], permutation[j] = permutation[j], permutation[i]
  end
  
  return permutation
end

local function UniformDeterministicSelection(parents, size)
  local perm_idx, permutation = 1, {}
  local offsprings = {}
  while #offsprings < size do
    if perm_idx > #permutation then
      permutation = GenPermutation(#parents)
      perm_idx = 1
    end
    local parent_idx = permutation[perm_idx]
    perm_idx = perm_idx + 1
    table.insert(offsprings, parents[parent_idx])
  end
  
  return offsprings
end

local function Crossover_1Point_2Offsprings(parent1, parent2)
  local xsite = math.random(1, #parent1)
  local offspring1, offspring2 = {}, {}
  for i = 1, xsite do
    offspring1[i], offspring2[i] = parent1[i], parent2[i]
  end
  for i = xsite + 1, #parent1 do
    offspring1[i], offspring2[i] = parent2[i], parent1[i]
  end
  
  return offspring1, offspring2
end

local function ProduceOffsprings2p(parent1, parent2, xsite1, xsite2)
  if xsite1 > xsite2 then
    xsite1, xsite2 = xsite2, xsite1
  end
  
  local offspring1, offspring2 = {}, {}
  for i = 1, xsite1 do
    offspring1[i], offspring2[i] = parent1[i], parent2[i]
  end
  for i = xsite1 + 1, xsite2 do
    offspring1[i], offspring2[i] = parent2[i], parent1[i]
  end
  for i = xsite2 + 1, #parent1 do
    offspring1[i], offspring2[i] = parent1[i], parent2[i]
  end
  
  return offspring1, offspring2
end

local function Crossover_12Point_2Offsprings(parent1, parent2)
  local xsite1 = math.random(1, #parent1 - 1)
  local xsite2 = math.random(xsite1 + 1, #parent1)
  
  return ProduceOffsprings2p(parent1, parent2, xsite1, xsite2)
end

local function Crossover_2Point_2Offsprings(parent1, parent2)
  local xsite1 = math.random(1, #parent1 - 2)
  local xsite2 = math.random(xsite1 + 1, #parent1 - 1)
  
  return ProduceOffsprings2p(parent1, parent2, xsite1, xsite2)
end

local function Crossover_Uniform_2Offsprings(parent1, parent2)
  local offspring1, offspring2 = {}, {}
  for i, allele in ipairs(parent1) do
    if math.random() < 0.5 then
      offspring1[i], offspring2[i] = allele, parent2[i]
    else
      offspring1[i], offspring2[i] = parent2[i], allele
    end
  end
  
  return offspring1, offspring2
end

local function Crossover_2Point_1Offspring(parent1, parent2)
  local xsite1 = math.random(1, #parent1 - 1)
  local xsite2 = math.random(xsite1 + 1, #parent1)
  local offspring = ProduceOffsprings2p(parent1, parent2, xsite1, xsite2)
  
  return offspring
end

local function Mutation(genotype, rate)
  for i, allele in ipairs(genotype) do
    if math.random() < rate then
      genotype[i] = (allele == 'A') and 'B' or 'A'
    end
  end
end

local function MutationCreator(rate)
  return function(genotype) return Mutation(genotype, rate) end
end

local function Crossover_Blend(parent1, parent2, r)
  local offspring = {}
  for i, allele in ipairs(parent1) do
    offspring[i] = r * allele + (1.0 - r) * parent2[i]
  end
  
  return offspring
end

local function CrossoverBlendCreator(r)
  return function(parent1, parent2) return Crossover_Blend(parent1, parent2, r) end
end

local function GausMutationChrom(genotype, genomes, step)
  local used = {}
  for k = 1, genomes do
    local loci
    while not loci or used[loci] do
      loci = math.random(1, #genotype)
    end
    used[loci] = true
    genotype[loci] = genotype[loci] + GausMutation(step)
  end
end

local function GausMutationCreator(genomes, step)
  return function(genotype) return GausMutationChrom(genotype, genomes, step) end
end

local function EV(m, n, init_pop, on_new_gen, crossover, mutation, selection_parent, selection_survival, overlapping)
  local parents = init_pop(m)
  if on_new_gen then
    on_new_gen(0, parents)
  end
  for gen = 1, MAX_GENERATIONS do
    local pool = selection_parent(parents, n)
    local offsprings = {}
    if crossover then
      for i = 1, #pool, 2 do
        local offspring1, offspring2 = crossover(pool[i], pool[i + 1])
        table.insert(offsprings, offspring1)
        if offspring2 then
          table.insert(offsprings, offspring2)
        end
      end
    else
      offsprings = pool
    end
    if mutation then
      for _, offspring in ipairs(offsprings) do
        mutation(offspring)
      end
    end
    if selection_survival then
      if overlapping then
        for _, genotype in ipairs(parents) do
          table.insert(offsprings, genotype)
        end
      end
      parents = selection_survival(offsprings, m)
    else
      parents = offsprings
    end
    if on_new_gen then
      on_new_gen(gen, parents)
    end
  end
  
  return gen, parents
end

local s_Seeds = {}
for run = 1, Max(RUNS, 50) do
  s_Seeds[run] = math.random(1, 100000)
end

local function reg_gen(points, gen, pop, run, rms_func)
  rms_func = rms_func or RMS
  
  local rms = rms_func(pop)
  local entry = points[gen] or {x = gen, y = 0}
  points[gen] = entry
  entry.y = entry.y + rms
end

local function CompareCrossoversDeterministic(filename)
  local graphs = {name_x = "Generations", name_y = "RMS Distance from Robbins Equilibrium", funcs = {}}
  local name_1p_2o = string.format("1-P,2 Offsprings")
  local name_12p_2o = string.format("1-2-P,2 Offsprings")
  local name_2p_2o = string.format("2-P,2 Offsprings")
  local name_2p_1o = string.format("2-P,1 Offspring")
  local name_uni_2o = string.format("Uniform(0.5),2 Offsprings")
  local points_1p_2o = {color = RGB_GREEN}
  local points_12p_2o = {color = RGB_CYAN}
  local points_2p_2o = {color = RGB_MAGENTA}
  local points_2p_1o = {color = RGB_ORANGE}
  local points_uni_2o = {color = RGB_RED}
  graphs.funcs[name_1p_2o] = points_1p_2o
  graphs.funcs[name_12p_2o] = points_12p_2o
  graphs.funcs[name_2p_2o] = points_2p_2o
  graphs.funcs[name_2p_1o] = points_2p_1o
  graphs.funcs[name_uni_2o] = points_uni_2o
  
  for run = 1, RUNS do
    math.randomseed(s_Seeds[run])
    EV(POPULATION_SIZE, 2 * POPULATION_SIZE, GenInitEdgePop, function(gen, pop)
      reg_gen(points_2p_1o, gen, pop, run)
    end, Crossover_2Point_1Offspring, nil, UniformDeterministicSelection)
  end
  
  for run = 1, RUNS do
    math.randomseed(s_Seeds[run])
    EV(POPULATION_SIZE, POPULATION_SIZE, GenInitEdgePop, function(gen, pop)
      reg_gen(points_1p_2o, gen, pop, run)
    end, Crossover_1Point_2Offsprings, nil, UniformDeterministicSelection)
  end
  
  for run = 1, RUNS do
    math.randomseed(s_Seeds[run])
    EV(POPULATION_SIZE, POPULATION_SIZE, GenInitEdgePop, function(gen, pop)
      reg_gen(points_12p_2o, gen, pop, run)
    end, Crossover_12Point_2Offsprings, nil, UniformDeterministicSelection)
  end
  
  for run = 1, RUNS do
    math.randomseed(s_Seeds[run])
    EV(POPULATION_SIZE, POPULATION_SIZE, GenInitEdgePop, function(gen, pop)
      reg_gen(points_2p_2o, gen, pop, run)
    end, Crossover_2Point_2Offsprings, nil, UniformDeterministicSelection)
  end
  
  for run = 1, RUNS do
    math.randomseed(s_Seeds[run])
    EV(POPULATION_SIZE, POPULATION_SIZE, GenInitEdgePop, function(gen, pop)
      reg_gen(points_uni_2o, gen, pop, run)
    end, Crossover_Uniform_2Offsprings, nil, UniformDeterministicSelection)
  end
  
  NormalizeGraphs(graphs, RUNS)
  
  local bmp = Bitmap.new(IMAGE_WIDTH, IMAGE_HEIGHT, RGB_BLACK)
  DrawGraphs(bmp, graphs, {int_x = true, skip_KP = true, min_x = 0, min_y = 0, div_x = 5, div_y = 4})
  bmp:WriteBMP(filename)
end

local function CompareCrossoversStochastic(filename)
  local graphs = {name_x = "Generations", name_y = "RMS Distance from Robbins Equilibrium", funcs = {}}
  local name_1p_2o = string.format("1-P,2 Offsprings")
  local name_12p_2o = string.format("1-2-P,2 Offsprings")
  local name_uni_2o = string.format("Uniform(0.5),2 Offsprings")
  local points_1p_2o = {color = RGB_GREEN}
  local points_12p_2o = {color = RGB_CYAN}
  local points_uni_2o = {color = RGB_RED}
  graphs.funcs[name_1p_2o] = points_1p_2o
  graphs.funcs[name_12p_2o] = points_12p_2o
  graphs.funcs[name_uni_2o] = points_uni_2o
  
  for run = 1, RUNS do
    math.randomseed(s_Seeds[run])
    EV(POPULATION_SIZE, POPULATION_SIZE, GenInitEdgePop, function(gen, pop)
      reg_gen(points_1p_2o, gen, pop, run)
    end, Crossover_1Point_2Offsprings, nil, UniformStochasticSelection)
  end
  
  for run = 1, RUNS do
    math.randomseed(s_Seeds[run])
    EV(POPULATION_SIZE, POPULATION_SIZE, GenInitEdgePop, function(gen, pop)
      reg_gen(points_12p_2o, gen, pop, run)
    end, Crossover_12Point_2Offsprings, nil, UniformStochasticSelection)
  end
  
  for run = 1, RUNS do
    math.randomseed(s_Seeds[run])
    EV(POPULATION_SIZE, POPULATION_SIZE, GenInitEdgePop, function(gen, pop)
      reg_gen(points_uni_2o, gen, pop, run)
    end, Crossover_Uniform_2Offsprings, nil, UniformStochasticSelection)
  end
  
  NormalizeGraphs(graphs, RUNS)
  
  local bmp = Bitmap.new(IMAGE_WIDTH, IMAGE_HEIGHT, RGB_BLACK)
  DrawGraphs(bmp, graphs, {int_x = true, skip_KP = true, min_x = 0, min_y = 0, div_x = 5, div_y = 4})
  bmp:WriteBMP(filename)
end

local function CompareMutations(filename)
  local graphs = {name_x = "Generations", name_y = "RMS Distance from Robbins Equilibrium", funcs = {}}
  local name_1L = string.format("1/L Mutation Rate")
  local name_2L = string.format("2/L Mutation Rate")
  local name_12L = string.format("1/2L Mutation Rate")
  local points_1L = {color = RGB_GREEN}
  local points_2L = {color = RGB_CYAN}
  local points_12L = {color = RGB_RED}
  graphs.funcs[name_1L] = points_1L
  graphs.funcs[name_2L] = points_2L
  graphs.funcs[name_12L] = points_12L
  
  for run = 1, RUNS do
    math.randomseed(s_Seeds[run])
    EV(POPULATION_SIZE, POPULATION_SIZE, GenInitUniformPop, function(gen, pop)
      reg_gen(points_1L, gen, pop, run)
    end, nil, MutationCreator(1.0 / GENOME_BITS), UniformDeterministicSelection)
  end
  
  for run = 1, RUNS do
    math.randomseed(s_Seeds[run])
    EV(POPULATION_SIZE, POPULATION_SIZE, GenInitUniformPop, function(gen, pop)
      reg_gen(points_2L, gen, pop, run)
    end, nil, MutationCreator(2.0 / GENOME_BITS), UniformDeterministicSelection)
  end
  
  for run = 1, RUNS do
    math.randomseed(s_Seeds[run])
    EV(POPULATION_SIZE, POPULATION_SIZE, GenInitUniformPop, function(gen, pop)
      reg_gen(points_12L, gen, pop, run)
    end, nil, MutationCreator(1.0 / (2 * GENOME_BITS)), UniformDeterministicSelection)
  end
  
  NormalizeGraphs(graphs, RUNS)
  
  local bmp = Bitmap.new(IMAGE_WIDTH, IMAGE_HEIGHT, RGB_BLACK)
  DrawGraphs(bmp, graphs, {int_x = true, skip_KP = true, min_x = 0, min_y = 0, div_x = 5, div_y = 4})
  bmp:WriteBMP(filename)
end

local function CompareMutationsHamming(RUNS, filename)
  local graphs = {name_x = "Generations", name_y = "Diversity", funcs = {}}
  local name_1L = string.format("1/L Mutation Rate")
  local name_2L = string.format("2/L Mutation Rate")
  local name_12L = string.format("1/2L Mutation Rate")
  local name_4L = string.format("4/L Mutation Rate")
  local points_1L = {color = RGB_GREEN}
  local points_2L = {color = RGB_CYAN}
  local points_12L = {color = RGB_RED}
  local points_4L = {color = RGB_ORANGE}
  graphs.funcs[name_1L] = points_1L
  graphs.funcs[name_2L] = points_2L
  graphs.funcs[name_12L] = points_12L
  graphs.funcs[name_4L] = points_4L
  
  for run = 1, RUNS do
    math.randomseed(s_Seeds[run])
    EV(POPULATION_SIZE, POPULATION_SIZE, GenInitRandomPop, function(gen, pop)
      reg_gen(points_1L, gen, pop, run, RMS_Hamming)
      if gen == MAX_GENERATIONS then
        print(string.format("Run: %d/%d, gen %d, RMS Hamming: %.2f", run, RUNS, gen, RMS_Hamming(pop)))
      end
    end, nil, MutationCreator(1.0 / GENOME_BITS), UniformDeterministicSelection)
  end
  
  for run = 1, RUNS do
    math.randomseed(s_Seeds[run])
    EV(POPULATION_SIZE, POPULATION_SIZE, GenInitRandomPop, function(gen, pop)
      reg_gen(points_2L, gen, pop, run, RMS_Hamming)
      if gen == MAX_GENERATIONS then
        print(string.format("Run: %d/%d, gen %d, RMS Hamming: %.2f", run, RUNS, gen, RMS_Hamming(pop)))
      end
    end, nil, MutationCreator(2.0 / GENOME_BITS), UniformDeterministicSelection)
  end
  
  for run = 1, RUNS do
    math.randomseed(s_Seeds[run])
    EV(POPULATION_SIZE, POPULATION_SIZE, GenInitRandomPop, function(gen, pop)
      reg_gen(points_12L, gen, pop, run, RMS_Hamming)
      if gen == MAX_GENERATIONS then
        print(string.format("Run: %d/%d, gen %d, RMS Hamming: %.2f", run, RUNS, gen, RMS_Hamming(pop)))
      end
    end, nil, MutationCreator(1.0 / (2 * GENOME_BITS)), UniformDeterministicSelection)
  end
  
  for run = 1, RUNS do
    math.randomseed(s_Seeds[run])
    EV(POPULATION_SIZE, POPULATION_SIZE, GenInitRandomPop, function(gen, pop)
      reg_gen(points_4L, gen, pop, run, RMS_Hamming)
      if gen == MAX_GENERATIONS then
        print(string.format("Run: %d/%d, gen %d, RMS Hamming: %.2f", run, RUNS, gen, RMS_Hamming(pop)))
      end
    end, nil, MutationCreator(4.0 / GENOME_BITS), UniformDeterministicSelection)
  end
  
  NormalizeGraphs(graphs, RUNS)
  
  local bmp = Bitmap.new(IMAGE_WIDTH, IMAGE_HEIGHT, RGB_BLACK)
  DrawGraphs(bmp, graphs, {int_x = true, skip_KP = true, min_x = 0, min_y = 10, div_x = 4, div_y = 7, max_y = 24})
  bmp:WriteBMP(filename)
end

local function CompareCrossoverMutastionsDiversity(filename)
  local graphs = {name_x = "Generations", name_y = "RMS Distance from Robbins Equilibrium", funcs = {}}
  local name_cm = string.format("Uniform Crossover and 1/L Mutation Rate")
  local name_c = string.format("Uniform Crossover Alone")
  local name_m = string.format("1/L Mutation Rate Alone")
  local points_cm = {color = RGB_GREEN}
  local points_c = {color = RGB_CYAN}
  local points_m = {color = RGB_RED}
  graphs.funcs[name_cm] = points_cm
  graphs.funcs[name_c] = points_c
  graphs.funcs[name_m] = points_m
  
  for run = 1, RUNS do
    math.randomseed(s_Seeds[run])
    EV(POPULATION_SIZE, POPULATION_SIZE, GenInitEdgePop, function(gen, pop)
      reg_gen(points_cm, gen, pop, run)
    end, Crossover_Uniform_2Offsprings, MutationCreator(1.0 / GENOME_BITS), UniformStochasticSelection)
  end
  
  for run = 1, RUNS do
    math.randomseed(s_Seeds[run])
    EV(POPULATION_SIZE, POPULATION_SIZE, GenInitEdgePop, function(gen, pop)
      reg_gen(points_c, gen, pop, run)
    end, Crossover_Uniform_2Offsprings, nil, UniformStochasticSelection)
  end
  
  for run = 1, RUNS do
    math.randomseed(s_Seeds[run])    
    EV(POPULATION_SIZE, POPULATION_SIZE, GenInitEdgePop, function(gen, pop)
      reg_gen(points_m, gen, pop, run)
    end, nil, MutationCreator(1 / GENOME_BITS), UniformStochasticSelection)
  end
  
  NormalizeGraphs(graphs, RUNS)
  
  local bmp = Bitmap.new(IMAGE_WIDTH, IMAGE_HEIGHT, RGB_BLACK)
  DrawGraphs(bmp, graphs, {int_x = true, skip_KP = true, min_x = 0, min_y = 0, div_x = 4, div_y = 7})
  bmp:WriteBMP(filename)
end

local function CompareDiversity(RUNS, filename)
  local graphs = {name_x = "Generations", name_y = "Diversity", funcs = {}}
  local name_c1L = string.format("Uniform Crossover,1/L Mutation Rate")
  local name_1L = string.format("1/L Mutation Rate")
  local name_c2L = string.format("Uniform Crossover,2/L Mutation Rate")
  local name_2L = string.format("2/L Mutation Rate")
  local points_c1L = {color = RGB_GREEN}
  local points_1L = {color = RGB_CYAN}
  local points_c2L = {color = RGB_RED}
  local points_2L = {color = RGB_ORANGE}
  graphs.funcs[name_c1L] = points_c1L
  graphs.funcs[name_1L] = points_1L
  graphs.funcs[name_c2L] = points_c2L
  graphs.funcs[name_2L] = points_2L
  
  for run = 1, RUNS do
    math.randomseed(s_Seeds[run])
    EV(POPULATION_SIZE, POPULATION_SIZE, GenInitRandomPop, function(gen, pop)
      reg_gen(points_c1L, gen, pop, run, RMS_Hamming)
      --if gen == MAX_GENERATIONS then
        print(string.format("Run: %d/%d, gen %d, RMS Hamming: %.2f", run, RUNS, gen, RMS_Hamming(pop)))
      --end
    end, Crossover_Uniform_2Offsprings, MutationCreator(1.0 / GENOME_BITS), UniformDeterministicSelection)
  end
  
  for run = 1, RUNS do
    math.randomseed(s_Seeds[run])
    EV(POPULATION_SIZE, POPULATION_SIZE, GenInitRandomPop, function(gen, pop)
      reg_gen(points_1L, gen, pop, run, RMS_Hamming)
      if gen == MAX_GENERATIONS then
        print(string.format("Run: %d/%d, gen %d, RMS Hamming: %.2f", run, RUNS, gen, RMS_Hamming(pop)))
      end
    end, nil, MutationCreator(1.0 / GENOME_BITS), UniformDeterministicSelection)
  end
  
  for run = 1, RUNS do
    math.randomseed(s_Seeds[run])
    EV(POPULATION_SIZE, POPULATION_SIZE, GenInitRandomPop, function(gen, pop)
      reg_gen(points_c2L, gen, pop, run, RMS_Hamming)
      if gen == MAX_GENERATIONS then
        print(string.format("Run: %d/%d, gen %d, RMS Hamming: %.2f", run, RUNS, gen, RMS_Hamming(pop)))
      end
    end, Crossover_Uniform_2Offsprings, MutationCreator(2.0 / GENOME_BITS), UniformDeterministicSelection)
  end
  
  for run = 1, RUNS do
    math.randomseed(s_Seeds[run])
    EV(POPULATION_SIZE, POPULATION_SIZE, GenInitRandomPop, function(gen, pop)
      reg_gen(points_2L, gen, pop, run, RMS_Hamming)
      if gen == MAX_GENERATIONS then
        print(string.format("Run: %d/%d, gen %d, RMS Hamming: %.2f", run, RUNS, gen, RMS_Hamming(pop)))
      end
    end, nil, MutationCreator(2.0 / GENOME_BITS), UniformDeterministicSelection)
  end
  
  NormalizeGraphs(graphs, RUNS)
  
  local bmp = Bitmap.new(IMAGE_WIDTH, IMAGE_HEIGHT, RGB_BLACK)
  DrawGraphs(bmp, graphs, {int_x = true, skip_KP = true, min_x = 0, min_y = 15, div_x = 4, div_y = 7, max_y = 24})
  bmp:WriteBMP(filename)
end

local function DispersionGausMutationDeterministic(filename)
  local graphs = {name_x = "Generations(Uniform-Deterministic Parent Selection)", name_y = "Population Dispersion", funcs = {}}
  local name_11 = string.format("1/10 GM step 1.0")
  local name_21 = string.format("2/10 GM step 1.0")
  local name_12 = string.format("1/10 GM step 2.0")
  local name_22 = string.format("2/10 GM step 2.0")
  local points_11 = {color = RGB_GREEN}
  local points_21 = {color = RGB_CYAN}
  local points_12 = {color = RGB_RED}
  local points_22 = {color = RGB_ORANGE}
  graphs.funcs[name_11] = points_11
  graphs.funcs[name_21] = points_21
  graphs.funcs[name_12] = points_12
  graphs.funcs[name_22] = points_22
  
  for run = 1, RUNS do
    math.randomseed(s_Seeds[run])
    EV(POPULATION_SIZE, POPULATION_SIZE, GenInitRealRandomPop, function(gen, pop)
      reg_gen(points_11, gen, pop, run, Dispersion)
    end, nil, GausMutationCreator(1, 1.0), UniformDeterministicSelection)
  end
  
  for run = 1, RUNS do
    math.randomseed(s_Seeds[run])
    EV(POPULATION_SIZE, POPULATION_SIZE, GenInitRealRandomPop, function(gen, pop)
      reg_gen(points_21, gen, pop, run, Dispersion)
    end, nil, GausMutationCreator(2, 1.0), UniformDeterministicSelection)
  end
  
  for run = 1, RUNS do
    math.randomseed(s_Seeds[run])
    EV(POPULATION_SIZE, POPULATION_SIZE, GenInitRealRandomPop, function(gen, pop)
      reg_gen(points_12, gen, pop, run, Dispersion)
    end, nil, GausMutationCreator(1, 2.0), UniformDeterministicSelection)
  end
  
  for run = 1, RUNS do
    math.randomseed(s_Seeds[run])
    EV(POPULATION_SIZE, POPULATION_SIZE, GenInitRealRandomPop, function(gen, pop)
      reg_gen(points_22, gen, pop, run, Dispersion)
    end, nil, GausMutationCreator(2, 2.0), UniformDeterministicSelection)
  end
  
  NormalizeGraphs(graphs, RUNS)
  
  local bmp = Bitmap.new(IMAGE_WIDTH, IMAGE_HEIGHT, RGB_BLACK)
  DrawGraphs(bmp, graphs, {int_x = true, skip_KP = true, min_x = 0, min_y = 0, div_x = 4, div_y = 6, max_y = 60})
  bmp:WriteBMP(filename)
end

local function DispersionGausMutationStochastic(filename)
  local graphs = {name_x = "Generations(Uniform-Stochastic Parent Selection)", name_y = "Population Dispersion", funcs = {}}
  local name_11 = string.format("1/10 GM step 1.0")
  local name_12 = string.format("1/10 GM step 2.0")
  local name_101 = string.format("10/10 GM step 1.0")
  local name_102 = string.format("10/10 GM step 2.0")
  local points_11 = {color = RGB_GREEN}
  local points_12 = {color = RGB_CYAN}
  local points_101 = {color = RGB_RED}
  local points_102 = {color = RGB_ORANGE}
  graphs.funcs[name_11] = points_11
  graphs.funcs[name_12] = points_12
  graphs.funcs[name_101] = points_101
  graphs.funcs[name_102] = points_102
  
  for run = 1, RUNS do
    math.randomseed(s_Seeds[run])
    EV(POPULATION_SIZE, POPULATION_SIZE, GenInitRealRandomPop, function(gen, pop)
      reg_gen(points_11, gen, pop, run, Dispersion)
    end, nil, GausMutationCreator(1, 1.0), UniformStochasticSelection)
  end
  
  for run = 1, RUNS do
    math.randomseed(s_Seeds[run])
    EV(POPULATION_SIZE, POPULATION_SIZE, GenInitRealRandomPop, function(gen, pop)
      reg_gen(points_12, gen, pop, run, Dispersion)
    end, nil, GausMutationCreator(1, 2.0), UniformStochasticSelection)
  end
  
  for run = 1, RUNS do
    math.randomseed(s_Seeds[run])
    EV(POPULATION_SIZE, POPULATION_SIZE, GenInitRealRandomPop, function(gen, pop)
      reg_gen(points_101, gen, pop, run, Dispersion)
    end, nil, GausMutationCreator(10, 1.0), UniformStochasticSelection)
  end
  
  for run = 1, RUNS do
    math.randomseed(s_Seeds[run])
    EV(POPULATION_SIZE, POPULATION_SIZE, GenInitRealRandomPop, function(gen, pop)
      reg_gen(points_102, gen, pop, run, Dispersion)
    end, nil, GausMutationCreator(10, 2.0), UniformStochasticSelection)
  end
  
  NormalizeGraphs(graphs, RUNS)
  
  local bmp = Bitmap.new(IMAGE_WIDTH, IMAGE_HEIGHT, RGB_BLACK)
  DrawGraphs(bmp, graphs, {int_x = true, skip_KP = true, min_x = 0, min_y = 0, div_x = 4, div_y = 6, max_y = 60})
  bmp:WriteBMP(filename)
end

local function DispersionBlend(filename)
  local graphs = {name_x = "Generations(Uniform-Deterministic Parent Selection)", name_y = "Population Dispersion", funcs = {}}
  local name_01 = string.format("Blend r=0.1")
  local name_03 = string.format("Blend r=0.3")
  local name_05 = string.format("Blend r=0.5")
  local points_01 = {color = RGB_GREEN}
  local points_03 = {color = RGB_CYAN}
  local points_05 = {color = RGB_RED}
  graphs.funcs[name_01] = points_01
  graphs.funcs[name_03] = points_03
  graphs.funcs[name_05] = points_05
  
  for run = 1, RUNS do
    math.randomseed(s_Seeds[run])
    EV(POPULATION_SIZE, POPULATION_SIZE, GenInitRealRandomPop, function(gen, pop)
      reg_gen(points_01, gen, pop, run, Dispersion)
    end, CrossoverBlendCreator(0.1), nil, UniformDeterministicSelection)
  end
  
  for run = 1, RUNS do
    math.randomseed(s_Seeds[run])
    EV(POPULATION_SIZE, POPULATION_SIZE, GenInitRealRandomPop, function(gen, pop)
      reg_gen(points_03, gen, pop, run, Dispersion)
    end, CrossoverBlendCreator(0.3), nil, UniformDeterministicSelection)
  end
  
  for run = 1, RUNS do
    math.randomseed(s_Seeds[run])
    EV(POPULATION_SIZE, POPULATION_SIZE, GenInitRealRandomPop, function(gen, pop)
      reg_gen(points_05, gen, pop, run, Dispersion)
    end, CrossoverBlendCreator(0.5), nil, UniformDeterministicSelection)
  end
  
  NormalizeGraphs(graphs, RUNS)
  
  local bmp = Bitmap.new(IMAGE_WIDTH, IMAGE_HEIGHT, RGB_BLACK)
  DrawGraphs(bmp, graphs, {int_x = true, skip_KP = true, min_x = 0, min_y = 0, div_x = 5, div_y = 5, max_y = 10})
  bmp:WriteBMP(filename)
end

local function DispersionBlendGausMutationDeterministic(filename)
  local graphs = {name_x = "Generations(Uniform-Deterministic Parent Selection)", name_y = "Population Dispersion", funcs = {}}
  local name_1 = string.format("Blend(0.5) 1/10 GM step 2.0")
  local name_2 = string.format("Blend(0.5) 2/10 GM step 2.0")
  local name_3 = string.format("Blend(0.5) 5/10 GM step 2.0")
  local name_4 = string.format("Blend(0.5) 10/10 GM step 2.0")
  local points_1 = {color = RGB_GREEN}
  local points_2 = {color = RGB_CYAN}
  local points_3 = {color = RGB_RED}
  local points_4 = {color = RGB_ORANGE}
  graphs.funcs[name_1] = points_1
  graphs.funcs[name_2] = points_2
  graphs.funcs[name_3] = points_3
  graphs.funcs[name_4] = points_4
  
  for run = 1, RUNS do
    math.randomseed(s_Seeds[run])
    EV(POPULATION_SIZE, POPULATION_SIZE, GenInitRealRandomPop, function(gen, pop)
      reg_gen(points_1, gen, pop, run, Dispersion)
    end, CrossoverBlendCreator(0.5), GausMutationCreator(1, 2.0), UniformDeterministicSelection)
  end
  
  for run = 1, RUNS do
    math.randomseed(s_Seeds[run])
    EV(POPULATION_SIZE, POPULATION_SIZE, GenInitRealRandomPop, function(gen, pop)
      reg_gen(points_2, gen, pop, run, Dispersion)
    end, CrossoverBlendCreator(0.5), GausMutationCreator(2, 2.0), UniformDeterministicSelection)
  end
  
  for run = 1, RUNS do
    math.randomseed(s_Seeds[run])
    EV(POPULATION_SIZE, POPULATION_SIZE, GenInitRealRandomPop, function(gen, pop)
      reg_gen(points_3, gen, pop, run, Dispersion)
    end, CrossoverBlendCreator(0.5), GausMutationCreator(5, 2.0), UniformDeterministicSelection)
  end
  
  for run = 1, RUNS do
    math.randomseed(s_Seeds[run])
    EV(POPULATION_SIZE, POPULATION_SIZE, GenInitRealRandomPop, function(gen, pop)
      reg_gen(points_4, gen, pop, run, Dispersion)
    end, CrossoverBlendCreator(0.5), GausMutationCreator(10, 2.0), UniformDeterministicSelection)
  end
  
  NormalizeGraphs(graphs, RUNS)
  
  local bmp = Bitmap.new(IMAGE_WIDTH, IMAGE_HEIGHT, RGB_BLACK)
  DrawGraphs(bmp, graphs, {int_x = true, skip_KP = true, min_x = 0, min_y = 0, div_x = 5, div_y = 8, max_y = 16})
  bmp:WriteBMP(filename)
end

local function CompareOverlappingVSNonOverlapping(filename)
  local graphs = {name_x = "Generations(Uniform-Deterministic Parent Selection)", name_y = "Population Dispersion", funcs = {}}
  local name_no_det= string.format("2-P,Non-Overlap,Det")
  local name_no_stoch = string.format("2-P,Non-Overlap,Stoch")
  local name_o_det= string.format("2-P,Overlap,Det")
  local name_o_stoch = string.format("2-P,Overlap,Stoch")
  local points_no_det = {color = RGB_GREEN}
  local points_no_stoch = {color = RGB_CYAN}
  local points_o_det = {color = RGB_RED}
  local points_o_stoch = {color = RGB_ORANGE}
  graphs.funcs[name_no_det] = points_no_det
  graphs.funcs[name_no_stoch] = points_no_stoch
  graphs.funcs[name_o_det] = points_o_det
  graphs.funcs[name_o_stoch] = points_o_stoch
  
  for run = 1, RUNS do
    math.randomseed(s_Seeds[run])
    EV(POPULATION_SIZE, POPULATION_SIZE, GenInitEdgePop, function(gen, pop)
      reg_gen(points_no_det, gen, pop, run, RMS)
    end, Crossover_12Point_2Offsprings, nil, UniformDeterministicSelection, UniformDeterministicSelection)
  end
  
  for run = 1, RUNS do
    math.randomseed(s_Seeds[run])
    EV(POPULATION_SIZE, POPULATION_SIZE, GenInitEdgePop, function(gen, pop)
      reg_gen(points_no_stoch, gen, pop, run, RMS)
    end, Crossover_12Point_2Offsprings, nil, UniformStochasticSelection, UniformDeterministicSelection)
  end
  
  for run = 1, RUNS do
    math.randomseed(s_Seeds[run])
    EV(POPULATION_SIZE, POPULATION_SIZE, GenInitEdgePop, function(gen, pop)
      reg_gen(points_o_det, gen, pop, run, RMS)
    end, Crossover_12Point_2Offsprings, nil, UniformDeterministicSelection, UniformStochasticSelection, "overlap")
  end
  
  for run = 1, RUNS do
    math.randomseed(s_Seeds[run])
    EV(POPULATION_SIZE, POPULATION_SIZE, GenInitEdgePop, function(gen, pop)
      reg_gen(points_o_stoch, gen, pop, run, RMS)
    end, Crossover_12Point_2Offsprings, nil, UniformStochasticSelection, UniformStochasticSelection, "overlap")
  end
  
  NormalizeGraphs(graphs, RUNS)
  
  local bmp = Bitmap.new(IMAGE_WIDTH, IMAGE_HEIGHT, RGB_BLACK)
  DrawGraphs(bmp, graphs, {int_x = true, skip_KP = true, min_x = 0, min_y = 0, div_x = 5, div_y = 4})
  bmp:WriteBMP(filename)
end

CompareCrossoversDeterministic(IMAGE_FILENAME_FIGURE_6_19)
CompareCrossoversStochastic(IMAGE_FILENAME_FIGURE_6_20)
CompareMutations(IMAGE_FILENAME_FIGURE_6_21)
MAX_GENERATIONS = 200 GENOME_BITS = 40
CompareMutationsHamming(50, IMAGE_FILENAME_FIGURE_6_22)
MAX_GENERATIONS = 100 GENOME_BITS = 5
CompareCrossoverMutastionsDiversity(IMAGE_FILENAME_FIGURE_6_23)
GENOME_BITS = 40
CompareDiversity(50, IMAGE_FILENAME_FIGURE_6_24)
GENOME_BITS = 5
MAX_GENERATIONS = 200
DispersionGausMutationDeterministic(IMAGE_FILENAME_FIGURE_6_25)
DispersionGausMutationStochastic(IMAGE_FILENAME_FIGURE_6_26)
MAX_GENERATIONS = 100
DispersionBlend(IMAGE_FILENAME_FIGURE_6_27)
DispersionBlendGausMutationDeterministic(IMAGE_FILENAME_FIGURE_6_28)
CompareOverlappingVSNonOverlapping(IMAGE_FILENAME_FIGURE_6_29)
