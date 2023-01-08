using CSV, DataFrames, DataFramesMeta

function filter_unchanged_cols(df)
    col_changed = map(eachcol(df)) do col
                      length(Set(col)) != 1
                  end
    return df[:,col_changed]
end

df = CSV.read("plots/df.csv", DataFrame)
df = filter_unchanged_cols(df) 
df_quench  = @subset(df, :Space_Time .== "quench")
df_rindler = @subset(df, :Space_Time .== "rindler")
df_flat    = @subset(df, :Space_Time .== "flat");

df_rindler


