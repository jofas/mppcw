using CSV, DataFrames, Statistics, Printf

function speedup_data()
  data = CSV.read("data.csv", delim=',', copycols=true)

  p1m = nothing

  for amnt_proc in groupby(data, :processes)
    p = amnt_proc[!, :processes][1]
    m = mean(amnt_proc[!, :time])

    if p == 1
      p1m = m
    end

    println(p, ",", p1m / m)
  end
end


function speedup_optimal()
  data = CSV.read("data.csv", delim=',', copycols=true)

  for amnt_proc in groupby(data, :processes)
    p = amnt_proc[!, :processes][1]
    println(p, ",", p)
  end
end


function avg_per_proc_and_size()
  data = CSV.read("data.csv", delim=',', copycols=true)

  gr = groupby(data, [:processes, :l])

  #l_line = string("\\diagbox{\$p\$}{\$n\$}", [string(" &", x[!, :l][1]) for x in
  #  [y for y in gr if y[!, :processes][1] == 1]]...)

  #println(l_line, "\\\\")
  #println("\\hline")

  cur_p = 0
  cur_line = ""

  for g in gr
    p, l = g[!, :processes][1], g[!, :l][1]

    m, s = mean(g[!, :time]), std(g[!, :time])

    if cur_p != p
      if cur_line != ""; println(cur_line, "\\\\"); end

      cur_p = p
      cur_line = "$p"
    end

    cur_line = cur_line * @sprintf " &%.3f &%.3f" m s
  end

  println(cur_line, "\\\\")
end


#speedup_data()
#speedup_optimal()
avg_per_proc_and_size()
