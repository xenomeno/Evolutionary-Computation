dofile("Statistics.lua")
dofile("Bitmap.lua")
dofile("Graphics.lua")
dofile("EC_Common.lua")
dofile("GA_Common.lua")

local RUNS                          = 5
local BIRTHS                        = 4000
local INIT_BOUND_MIN                = -10.0
local INIT_BOUND_MAX                = 10.0
local BOUND_MIN                     = -10.0
local BOUND_MAX                     = 10.0
local BINARY_PRECISION              = 0.00001
local GENOME_BITS                   = math.ceil(math.log((BOUND_MAX - BOUND_MIN) / BINARY_PRECISION, 2))
local BINARY_MAX                    = math.pow(2, GENOME_BITS) - 1
local CHROMOSOME_LENGTH             = 5 * GENOME_BITS
local MUTATION_STEP                 = 0.1
local PROB_CROSSOVER                = 1.0
local PROB_UNIFORM_CROSSOVER_RATE   = 0.2
local PROB_MUTATION                 = 1 / CHROMOSOME_LENGTH

local INTEGER_GENOMES               = 10
local INTEGER_GENOME_MAX            = 31
local INTEGER_TARGET_ALLELE         = 4
local INTEGER_MUTATION_STEP         = 1
local INTEGER_GENOME_BITS           = math.ceil(math.log(INTEGER_GENOME_MAX, 2))
local FITNESS_FUNCTION              = false

local GRAPH_NODES                   = 30
local GRAPH_EDGES                   = GRAPH_NODES * GRAPH_NODES
local GRAPH_BITS                    = math.ceil(math.log(GRAPH_EDGES, 2))

local IMAGE_WIDTH                   = 1000
local IMAGE_HEIGHT                  = 1000
local IMAGE_FILENAME_FIGURE_5_1     = "EC/Figure5_1.bmp"
local IMAGE_FILENAME_FIGURE_5_2     = "EC/Figure5_2.bmp"
local IMAGE_FILENAME_FIGURE_5_3     = "EC/Figure5_3.bmp"
local IMAGE_FILENAME_FIGURE_5_4     = "EC/Figure5_4.bmp"
local IMAGE_FILENAME_FIGURE_5_5     = "EC/Figure5_5.bmp"
local IMAGE_FILENAME_FIGURE_5_6     = "EC/Figure5_6.bmp"
local IMAGE_FILENAME_FIGURE_5_7     = "EC/Figure5_7.bmp"
local IMAGE_FILENAME_FIGURE_5_8     = "EC/Figure5_8.bmp"
local IMAGE_FILENAME_FIGURE_5_9     = "EC/Figure5_9.bmp"
local IMAGE_FILENAME_FIGURE_5_10    = "EC/Figure5_10.bmp"

local function RealToDiscrete(x)
  return math.floor(BINARY_MAX * (x - BOUND_MIN) / (BOUND_MAX - BOUND_MIN))
end

local function DescreteToReal(x)
  return BOUND_MIN + (BOUND_MAX - BOUND_MIN) * x / BINARY_MAX
end

local function IntegerToBitsring(x)
  local bitstring = {}
  while #bitstring < GENOME_BITS do
    table.insert(bitstring, ((x & 1) == 1) and "1" or "0")
    x = x >> 1
  end
  
  return string.reverse(table.concat(bitstring, ""))
end

local function BitstringToInteger(bitstring)
  local x, pow2 = 0, 1
  for i = #bitstring, 1, -1 do
    if string.sub(bitstring, i, i) == "1" then
      x = x + pow2
    end
    pow2 = pow2 << 1
  end
  
  return x
end

local function BinaryEncode(chrom)
  local d1 = RealToDiscrete(chrom.x1)
  local d2 = RealToDiscrete(chrom.x2)
  local d3 = RealToDiscrete(chrom.x3)
  local d4 = RealToDiscrete(chrom.x4)
  local d5 = RealToDiscrete(chrom.x5)
  
  local bitstring1 = IntegerToBitsring(d1)
  local bitstring2 = IntegerToBitsring(d2)
  local bitstring3 = IntegerToBitsring(d3)
  local bitstring4 = IntegerToBitsring(d4)
  local bitstring5 = IntegerToBitsring(d5)
  
  return PackBitstring(string.format("%s%s%s%s%s", bitstring1, bitstring2, bitstring3, bitstring4, bitstring5))
end

local function BinaryDecode(chrom)
  local bitstring = UnpackBitstring(chrom)
  
  local bitstring1 = string.sub(bitstring, 1, GENOME_BITS)
  local bitstring2 = string.sub(bitstring, GENOME_BITS + 1, 2 * GENOME_BITS)
  local bitstring3 = string.sub(bitstring, 2 * GENOME_BITS + 1, 3 * GENOME_BITS)
  local bitstring4 = string.sub(bitstring, 3 * GENOME_BITS + 1, 4 * GENOME_BITS)
  local bitstring5 = string.sub(bitstring, 4 * GENOME_BITS + 1, 5 * GENOME_BITS)
  
  local d1 = BitstringToInteger(bitstring1)
  local d2 = BitstringToInteger(bitstring2)
  local d3 = BitstringToInteger(bitstring3)
  local d4 = BitstringToInteger(bitstring4)
  local d5 = BitstringToInteger(bitstring5)
  
  local x1 = DescreteToReal(d1)
  local x2 = DescreteToReal(d2)
  local x3 = DescreteToReal(d3)
  local x4 = DescreteToReal(d4)
  local x5 = DescreteToReal(d5)
  
  return x1, x2, x3, x4, x5, bitstring1, bitstring2, bitstring3, bitstring4, bitstring5
end

local function BBF1(chrom, binary)
  local x1, x2, x3, x4, x5
  if binary then
    x1, x2, x3, x4, x5 = BinaryDecode(chrom)
  else
    x1, x2, x3, x4, x5 = chrom.x1, chrom.x2, chrom.x3, chrom.x4, chrom.x5
  end
  
  local v1 = 50 - x1 * x1
  local v2 = 50 - x2 * x2
  local v3 = 50 - x3 * x3
  local v4 = 50 - x4 * x4
  local v5 = 50 - x5 * x5

  return (v1 + v2 + v3 + v4 + v5) / 5
end

local function BBF2(chrom, binary)
  local x1, x2, x3, x4, x5
  if binary then
    x1, x2, x3, x4, x5 = BinaryDecode(chrom)
  else
    x1, x2, x3, x4, x5 = chrom.x1, chrom.x2, chrom.x3, chrom.x4, chrom.x5
  end
  
  return x1 * x1 + x2 * x2 + x3 * x3 + x4 * x4 + x5 * x5
end

local BBF3_Target = {}
for i = 1, INTEGER_GENOMES do
  BBF3_Target[i] = INTEGER_TARGET_ALLELE
end

local function BitstringToIntegerGenome(bitstring)
  local genome = {}
  for i = 1, INTEGER_GENOMES do
    local bit_pos = (i - 1) * INTEGER_GENOME_BITS + 1
    local bitstr = string.sub(bitstring, bit_pos, bit_pos + INTEGER_GENOME_BITS - 1)
    genome[i] = BitstringToInteger(bitstr)
  end
  
  return genome
end

local function IntegerGenomeEncode(genome)
  local bitstrings = {}
  for i, allele in ipairs(genome) do
    bitstrings[i] = IntegerToBitsring(allele)
  end
  local bitstring = table.concat(bitstrings, "")
  
  return PackBitstring(bitstring)
end

local function IntegerGenomeDecode(genome)
  local bitstring = UnpackBitstring(genome)
  return BitstringToIntegerGenome(bitstring)
end

local function BBF3(chrom, chrom_type)
  local genome = (chrom_type == "binary") and IntegerGenomeDecode(chrom) or chrom
  local dist2 = 0
  for i = 1, INTEGER_GENOMES do
    local delta = BBF3_Target[i] - genome[i]
    dist2 = dist2 + delta * delta
  end
  
  return -math.sqrt(dist2)
end

local function BBF4(chrom, chrom_type)
  local genome = (chrom_type == "binary") and IntegerGenomeDecode(chrom) or chrom
  local mismatched = 0
  for i = 1, INTEGER_GENOMES do
    mismatched = mismatched + ((BBF3_Target[i] ~= genome[i]) and 1 or 0)
  end
  
  return -mismatched
end

local BBF5_Target = {}
for row = 1, GRAPH_NODES do
  BBF5_Target[row] = {}
  local row_value = ((row & 1) == 1) and "1" or "0"
  for col = 1, GRAPH_NODES do
    BBF5_Target[row][col] = row_value
  end
end

local function Linearize(matrix, chrom_type)
  local bitstring
  if chrom_type == "row" then
    local by_rows = {}
    for row = 1, GRAPH_NODES do
      by_rows[row] = table.concat(matrix[row], "")
    end
    bitstring = table.concat(by_rows, "")
  else
    local by_cols = {}
    for col = 1, GRAPH_NODES do
      by_cols[col] = {}
      for row = 1, GRAPH_NODES do
        by_cols[col][row] = matrix[row][col]
      end
      by_cols[col] = table.concat(by_cols[col], "")
    end
    bitstring = table.concat(by_cols, "")
  end

  return PackBitstring(bitstring)
end

local function Unlinearize(chrom, chrom_type)
  local bitstring = UnpackBitstring(chrom)
  local matrix = {}
  local row, col = 1, 1
  if chrom_type == "row" then
    for i = 1, #bitstring do
      local value = string.sub(bitstring, i, i)
      matrix[row] = matrix[row] or {}
      matrix[row][col] = value
      col = col + 1
      if col > GRAPH_NODES then
        row = row + 1
        col = 1
      end
    end
  else
    for i = 1, #bitstring do
      local value = string.sub(bitstring, i, i)
      matrix[row] = matrix[row] or {}
      matrix[row][col] = value
      row = row + 1
      if row > GRAPH_NODES then
        col = col + 1
        row = 1
      end
    end
  end
  
  return matrix
end

local function BBF5(chrom, linearization)
  local matrix = Unlinearize(chrom, linearization)
  local match = 0
  for row = 1, GRAPH_NODES do
    for col = 1, GRAPH_NODES do
      match = match + ((matrix[row][col] == BBF5_Target[row][col]) and 1 or 0)
    end
  end
  
  return match
end

local function Copy(chrom)
  return {x1 = chrom.x1, x2 = chrom.x2, x3 = chrom.x3, x4 = chrom.x4, x5 = chrom.x5}
end

local function Mutation(chrom)
  local x1, x2, x3, x4, x5 = chrom.x2, chrom.x2, chrom.x3, chrom.x4, chrom.x5
  if math.random() < 0.5 then
    local new = x1 + GausMutation(MUTATION_STEP)
    if BOUND_MIN <= new and new <= BOUND_MAX then
      x1 = new
    end
  end
  if math.random() < 0.5 then
    local new = x2 + GausMutation(MUTATION_STEP)
    if BOUND_MIN <= new and new <= BOUND_MAX then
      x2 = new
    end
  end
  if math.random() < 0.5 then
    local new = x3 + GausMutation(MUTATION_STEP)
    if BOUND_MIN <= new and new <= BOUND_MAX then
      x3 = new
    end
  end
  if math.random() < 0.5 then
    local new = x4 + GausMutation(MUTATION_STEP)
    if BOUND_MIN <= new and new <= BOUND_MAX then
      x4 = new
    end
  end
  if math.random() < 0.5 then
    local new = x5 + GausMutation(MUTATION_STEP)
    if BOUND_MIN <= new and new <= BOUND_MAX then
      x5 = new
    end
  end
  
  return {x1 = x1, x2 = x2, x3 = x4, x4 = x4, x5 = x5}
end

local function EvalPop(pop)
  local min_fitness, max_fitness
  local total_fitness = 0.0
  local best_ind
  for _, ind in ipairs(pop) do
    local fitness = FITNESS_FUNCTION(ind.chrom, pop.chrom_type)
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

local function GenInitPop(size, on_birth, chrom_type)
  local pop = {chrom_type = chrom_type}
  while #pop < size do
    local chrom = {}
    if FITNESS_FUNCTION == BBF3 or FITNESS_FUNCTION == BBF4 then
      for i = 1, INTEGER_GENOMES do
        chrom[i] = math.random(0, INTEGER_GENOME_MAX)
      end
      if chrom_type == "binary" then
        chrom = IntegerGenomeEncode(chrom)
      end
    else
      chrom = {
          x1 = GenRand(INIT_BOUND_MIN, INIT_BOUND_MAX),
          x2 = GenRand(INIT_BOUND_MIN, INIT_BOUND_MAX),
          x3 = GenRand(INIT_BOUND_MIN, INIT_BOUND_MAX),
          x4 = GenRand(INIT_BOUND_MIN, INIT_BOUND_MAX),
          x5 = GenRand(INIT_BOUND_MIN, INIT_BOUND_MAX)
      }
      if chrom_type == "binary" then
        chrom = BinaryEncode(chrom)
      end
    end
    table.insert(pop, {chrom = chrom, birthdate = #pop + 1})
    if on_birth then
      EvalPop(pop)
      on_birth(pop, #pop)
    end
  end

  return pop
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

local function TruncateSelection(parents, offsprings, size)
  local pool = GetMergedSorted(parents, offsprings)
  while #pool > size do
    table.remove(pool)
  end
  
  return pool
end

local function EV(m, n, init_pop, on_birth, selection, no_mutation, max_births, non_overlapping)
  max_births = max_births or BIRTHS
  
  local parents = init_pop(m, on_birth)
  local birth = #parents
  local gen = 1
  local offsprings = {}
  while birth < max_births do
    local offsprings = {}
    for i = 1, n do
      local parent = parents[math.random(1, #parents)]
      local chrom = no_mutation and Copy(parent.chrom) or Mutation(parent.chrom)
      local offspring = {chrom = chrom}
      offspring.fitness = FITNESS_FUNCTION(offspring.chrom)
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

local function ES(m, n, init_pop_factory, on_birth)
  EV(m, n, init_pop_factory, on_birth, TruncateSelection)
end

local function BinaryTwoPointsCrossover(parent1, parent2)
  local bitstring1 = UnpackBitstring(parent1)
  local bitstring2 = UnpackBitstring(parent2)
  local xsite1 = math.random(1, CHROMOSOME_LENGTH - 1)
  local xsite2 = math.random(xsite1 + 1, CHROMOSOME_LENGTH)
  local offspring1, offspring2 = {}, {}
  for i = 1, CHROMOSOME_LENGTH do
    if ((i < xsite1) or (i >= xsite2)) then
      offspring1[i] = string.sub(bitstring1, i, i)
      offspring2[i] = string.sub(bitstring2, i, i)
    else
      offspring1[i] = string.sub(bitstring2, i, i)
      offspring2[i] = string.sub(bitstring1, i, i)
    end
  end
  
  return PackBitstring(table.concat(offspring1, "")), PackBitstring(table.concat(offspring2, ""))
end

local function BinaryOnePointCrossover(parent1, parent2)
  local offspring1, offspring2 = CopyBitstring(parent1), CopyBitstring(parent2)
  if math.random() < PROB_CROSSOVER then
    local xsite = math.random(1, parent1.bits)
    ExchangeTailBits(offspring1, offspring2, xsite)
  end

  return offspring1, offspring2
end

local function BinaryUniformCrossover(parent1, parent2)
  local bitstring1, bitstring2 = UnpackBitstring(parent1), UnpackBitstring(parent2)
  local offspring1, offspring2 = {}, {}
  local len = #bitstring1
  for i = 1, len do
    local allele1, allele2 = string.sub(bitstring1, i, i), string.sub(bitstring2, i, i)
    if math.random() < PROB_UNIFORM_CROSSOVER_RATE then
      offspring1[i], offspring2[i] = allele1, allele2
    else
      offspring1[i], offspring2[i] = allele2, allele1
    end
  end
  offspring1, offspring2 = table.concat(offspring1, ""), table.concat(offspring2, "")
  
  return PackBitstring(offspring1), PackBitstring(offspring2)
end

local function BinaryMutation(chrom)
  local word_idx, bit_pos, power2 = 1, 1, 1
  for bit = 1, chrom.bits do
    if FlipCoin(PROB_MUTATION) then
      local word = chrom[word_idx]
      local allele = word & power2
      chrom[word_idx] = (allele ~= 0) and (word - power2) or (word + power2)
    end
    bit_pos = bit_pos + 1
    power2 = power2 << 1
    if bit_pos > GetBitstringWordSize() then
      word_idx = word_idx + 1
      bit_pos, power2 = 1, 1
    end
  end
end

local function IntegerUniformCrossover(parent1, parent2)
  local offspring1, offspring2 = {}, {}
  for i, allele in ipairs(parent1) do
    if math.random() < PROB_UNIFORM_CROSSOVER_RATE then
      offspring1[i], offspring2[i] = allele, parent2[i]
    else
      offspring1[i], offspring2[i] = parent2[i], allele
    end
  end
  
  return offspring1, offspring2
end

local function IntegerMutation(chrom)
  for i, allele in ipairs(chrom) do
    local new = allele + GausMutation(INTEGER_MUTATION_STEP)
    if BOUND_MIN <= new and new <= BOUND_MAX then
      chrom[i] = math.floor(new + 0.5)
    end
  end
end

local function SymbolicMutation(chrom)
  for i, allele in ipairs(chrom) do
    chrom[i] = math.random(1, INTEGER_GENOMES)
  end
end

local function BinaryEV(m, n, init_pop, chrom_type, crossovers, mutations, on_birth)
  max_births = max_births or BIRTHS
  
  local crossover = crossovers[chrom_type]
  local mutation = mutations[chrom_type]
  
  local parents = init_pop(m, on_birth, chrom_type)
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
      mutation(offspring1)
      table.insert(offsprings, {chrom = offspring1, fitness = FITNESS_FUNCTION(offspring1, chrom_type)})
      if #offsprings < n then
        mutation(offspring2)
        table.insert(offsprings, {chrom = offspring2, fitness = FITNESS_FUNCTION(offspring2, chrom_type)})
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

local function GA(m, n, init_pop, chrom_type, crossovers, mutations, on_birth)
  BinaryEV(m, n, init_pop, chrom_type, crossovers, mutations, on_birth)
end

local s_Seeds = {}
for run = 1, RUNS do
  s_Seeds[run] = math.random(1, 100000)
end

local function FigureES(m, n, init_pop_factory, mutation_step, filename)
  MUTATION_STEP = mutation_step
  
  local names, graphs = {}, {name_x = "Number of Births", name_y = string.format("Optimum: Best-So-Far"), funcs = {}}
  local colors = {RGB_CYAN, RGB_GREEN, RGB_RED, RGB_MAGENTA, RGB_WHITE}
  for run = 1, RUNS do
    local name = string.format("ES(%d, %d) run %d", m, n, run)
    table.insert(names, name)
    graphs.funcs[name] = {color = colors[run]}
  end
  
  local best_so_far = {}
  print(string.format("Run\tOptimum\t\tParameter values"))
  for run = 1, RUNS do
    math.randomseed(s_Seeds[run])
    local best_ind
    ES(m, n, init_pop_factory, function(pop, birth)
      local name = names[run]
      best_ind = (not best_ind or pop.best_ind.fitness > best_ind.fitness) and pop.best_ind or best_ind
      points = graphs.funcs[name]
      if points[birth] then
        points[birth].y = points[birth].y + best_ind.fitness
      else
        points[birth] = {x = birth, y = best_ind.fitness}
      end
    end)
    best_so_far[run] = best_ind
    local chrom = best_ind.chrom
    local x1, x2, x3, x4, x5 = chrom.x1, chrom.x2, chrom.x3, chrom.x4, chrom.x5
    print(string.format("%2d\t%f\t%f %f %f %f %f", run, best_ind.fitness, x1, x2, x3, x4, x5))
  end
  
  for _, points in pairs(graphs.funcs) do
    for _, pt in ipairs(points) do
        pt.y = (pt.y < -50) and -50 or pt.y
    end
  end
  
  local bmp = Bitmap.new(IMAGE_WIDTH, IMAGE_HEIGHT, RGB_BLACK)
  DrawGraphs(bmp, graphs, {int_x = true, skip_KP = true, min_y = -50.0 })
  bmp:WriteBMP(filename)
end

local function FigureGA(m, n, init_pop_factory, crossovers, mutations, filename)
  local names, graphs = {}, {name_x = "Number of Births", name_y = string.format("Optimum: Best-So-Far"), funcs = {}}
  local colors = {RGB_CYAN, RGB_GREEN, RGB_RED, RGB_MAGENTA, RGB_WHITE}
  for run = 1, RUNS do
    local name = string.format("GA(%d, %d) run %d", m, n, run)
    table.insert(names, name)
    graphs.funcs[name] = {color = colors[run]}
  end
  
  local best_so_far = {}
  print(string.format("Run\tOptimum\t\tParameter values\t\t\t\t\t\tBit Strings"))
  for run = 1, RUNS do
    math.randomseed(s_Seeds[run])
    local best_ind
    GA(m, n, init_pop_factory, "binary", crossovers, mutations, function(pop, birth)
      local name = names[run]
      best_ind = (not best_ind or pop.best_ind.fitness > best_ind.fitness) and pop.best_ind or best_ind
      points = graphs.funcs[name]
      if points[birth] then
        points[birth].y = points[birth].y + best_ind.fitness
      else
        points[birth] = {x = birth, y = best_ind.fitness}
      end
    end)
    best_so_far[run] = best_ind
    if chrom_type then
      local genome = (chrom_type == "binary") and BitstringToIntegerGenome(best_ind.chrom) or best_ind.chrom
      print(string.format("%2d\t%f\t{%s}", run, best_ind.fitness, table.concat(genome, ", ")))
    else
      local x1, x2, x3, x4, x5, bs1, bs2, bs3, bs4, bs5 = BinaryDecode(best_ind.chrom)
      print(string.format("%2d\t%f\t%f %f %f %f %f\t\t\t%s %s %s %s %s", run, best_ind.fitness, x1, x2, x3, x4, x5, bs1, bs2, bs3, bs4, bs5))
    end
  end
  
  for _, points in pairs(graphs.funcs) do
    for _, pt in ipairs(points) do
        pt.y = (pt.y < -50) and -50 or pt.y
    end
  end
  
  local bmp = Bitmap.new(IMAGE_WIDTH, IMAGE_HEIGHT, RGB_BLACK)
  DrawGraphs(bmp, graphs, {int_x = true, skip_KP = true, min_y = -50.0})
  bmp:WriteBMP(filename)
end

local function FigureGAs(m, n, init_pop_factory, chrom_types, crossovers, mutations, filename)
  local names, graphs = {}, {name_x = "Number of Births", name_y = string.format("Optimum: Best-So-Far"), funcs = {}}
  local colors = {RGB_CYAN, RGB_GREEN, RGB_RED, RGB_MAGENTA, RGB_WHITE}
  for i, chrom_type in ipairs(chrom_types) do
    local name = string.format("GA(%d, %d) with %s genes", m, n, chrom_type)
    names[chrom_type] = name
    graphs.funcs[name] = {color = colors[i]}
  end
  
  for _, chrom_type in ipairs(chrom_types) do
    math.randomseed(s_Seeds[1])
    local best_ind
    GA(m, n, init_pop_factory, chrom_type, crossovers, mutations, function(pop, birth)
      local name = names[chrom_type]
      best_ind = (not best_ind or pop.best_ind.fitness > best_ind.fitness) and pop.best_ind or best_ind
      points = graphs.funcs[name]
      if points[birth] then
        points[birth].y = points[birth].y + best_ind.fitness
      else
        points[birth] = {x = birth, y = best_ind.fitness}
      end
    end)
    if chrom_type == "row" or chrom_type == "col" then
      local matrix = Unlinearize(best_ind.chrom, chrom_type)
      print(string.format("Graph (%dx%d) by %s %f(%.2f%%)\t\t\tTarget", GRAPH_NODES, GRAPH_NODES, chrom_type, best_ind.fitness, 100.0 * best_ind.fitness / GRAPH_EDGES))
      for row = 1, GRAPH_NODES do
        print(string.format("%s\t\t\t%s", table.concat(matrix[row], ""), string.rep(((row & 1) == 1) and "1" or "0", GRAPH_NODES)))
      end
    else
      local genome = (chrom_type == "binary") and IntegerGenomeDecode(best_ind.chrom) or best_ind.chrom
      print(string.format("%s\t%f\t{%s}", chrom_type, best_ind.fitness, table.concat(genome, ", ")))
    end
  end
  
  local min_y = (FITNESS_FUNCTION == BBF5) and 100 or -15
  local y_line = (FITNESS_FUNCTION ~= BBF5) and 0 or nil
  local div_y = (FITNESS_FUNCTION == BBF5) and 6 or 4
  
  for _, points in pairs(graphs.funcs) do
    for _, pt in ipairs(points) do
        pt.y = (pt.y < min_y) and min_y or pt.y
    end
  end
  
  local bmp = Bitmap.new(IMAGE_WIDTH, IMAGE_HEIGHT, RGB_BLACK)
  DrawGraphs(bmp, graphs, {int_x = true, skip_KP = true, min_y = min_y, y_line = y_line, div_y = div_y, int_y = true})
  bmp:WriteBMP(filename)
end

local function VectorInitPopFactory(init_min_bound, init_max_bound, min_bound, max_bound)
  return function(size, on_birth, chrom_type)
    INIT_BOUND_MIN = init_min_bound
    INIT_BOUND_MAX = init_max_bound
    BOUND_MIN = min_bound
    BOUND_MAX = max_bound
    
    return GenInitPop(size, on_birth, chrom_type)
  end
end

local function GenRandomMatrix()
  local matrix = {}
  for i = 1, GRAPH_NODES do
    matrix[i] = {}
    for j = 1, GRAPH_NODES do
      matrix[i][j] = (math.random() < 0.5) and "1" or "0"
    end
  end
  
  return matrix
end

local function GraphInitPopFactory()
  return function(size, on_birth, chrom_type)
    local pop = {chrom_type = chrom_type}
    while #pop < size do
      -- NOTE: no need to generate a matrix and linearize it since its random - just generating a random string will do the same
      --local graph = PackBitstring(GenRandomBitstring(GRAPH_EDGES))
      local matrix_adj = GenRandomMatrix()
      local graph = Linearize(matrix_adj, chrom_type)
      table.insert(pop, {chrom = graph, birthdate = #pop + 1})
      if on_birth then
        EvalPop(pop)
        on_birth(pop, #pop)
      end
    end

    return pop
  end
end

FITNESS_FUNCTION = BBF1
local crossovers = {["binary"] = BinaryTwoPointsCrossover, ["integer"] = IntegerUniformCrossover, ["symbolic"] = IntegerUniformCrossover}
local mutations = {["binary"] = BinaryMutation, ["integer"] = IntegerMutation, ["symbolic"] = SymbolicMutation}
FigureES(5, 25, VectorInitPopFactory(-10.0, 10.0, -10.0, 10.0), 0.1, IMAGE_FILENAME_FIGURE_5_1)
FigureES(15, 15, VectorInitPopFactory(-100.0, 100.0, -100.0, 100.0), 1.0, IMAGE_FILENAME_FIGURE_5_2)
FigureGA(15, 15, VectorInitPopFactory(-10.0, 10.0, -10.0, 10.0), crossovers, mutations, IMAGE_FILENAME_FIGURE_5_3)

FITNESS_FUNCTION = BBF2
FigureES(5, 25, VectorInitPopFactory(-4.0, 5.0, -4.0, 5.0), 0.1, IMAGE_FILENAME_FIGURE_5_4)
FigureES(15, 15, VectorInitPopFactory(-4.0, 5.0, -4.0, 5.0), 0.1, IMAGE_FILENAME_FIGURE_5_5)
FigureGA(15, 15, VectorInitPopFactory(-4.0, 5.0, -4.0, 5.0), crossovers, mutations, IMAGE_FILENAME_FIGURE_5_6)

FITNESS_FUNCTION = BBF3
BIRTHS = 5000
FigureGAs(50, 50, VectorInitPopFactory(0, 31, 0, 31), {"integer", "binary"}, crossovers, mutations, IMAGE_FILENAME_FIGURE_5_7)

FITNESS_FUNCTION = BBF4
FigureGAs(50, 50, VectorInitPopFactory(0, 31, 0, 31), {"symbolic", "binary"}, crossovers, mutations, IMAGE_FILENAME_FIGURE_5_8)

FITNESS_FUNCTION = BBF5
BIRTHS = 10000
PROB_MUTATION = 1.0 / GRAPH_EDGES
local crossovers = {["row"] = BinaryOnePointCrossover, ["col"] = BinaryOnePointCrossover}
local mutations = {["row"] = BinaryMutation, ["col"] = BinaryMutation}
FigureGAs(100, 100, GraphInitPopFactory(), {"row", "col"}, crossovers, mutations, IMAGE_FILENAME_FIGURE_5_9)
local crossovers = {["row"] = BinaryUniformCrossover, ["col"] = BinaryUniformCrossover}
FigureGAs(100, 100, GraphInitPopFactory(), {"row", "col"}, crossovers, mutations, IMAGE_FILENAME_FIGURE_5_10)