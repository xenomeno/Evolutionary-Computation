function GenRand(min, max)
  return min + math.random() * (max - min)
end

local s_GausNorm = 1.0 / math.sqrt(2.0 / math.pi)

function GausMutation(s)
  return GenGausNoise(0.0, s * s_GausNorm)
end

function PlotPop(pop, birth, points, param)
  if points[birth] then
    points[birth].y = points[birth].y + pop[param]
  else
    points[birth] = {x = birth, y = pop[param]}
  end
end

function NormalizeGraphs(graphs, runs, min_y)
  local scale = 1.0 / runs
  for name, func in pairs(graphs.funcs) do
    local start_idx = func[0] and 0 or 1
    local stddev_interval = func.stddev_interval
    local stddev_start = func.stddev_start
    for k = start_idx, #func do
      func[k].y = func[k].y * scale
      local values = func[k].values
      if stddev_interval and values and (k - stddev_start) % stddev_interval == 0 then
        local mean = 0.0
        for _, v in ipairs(values) do
          mean = mean + v
        end
        mean = mean / #values
        local variance = 0.0
        for _, v in ipairs(values) do
          variance = variance + (v - mean) * (v - mean)
        end
        variance = variance / #values
        func[k].stddev = math.sqrt(variance)
      end
      if min_y then
        func[k].y = (func[k].y > min_y) and func[k].y or min_y
      end
    end
  end
end

