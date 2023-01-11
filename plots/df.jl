using CSV, DataFrames, DataFramesMeta

function filter_unchanged_cols(df)
    col_changed = map(eachcol(df)) do col
                      length(Set(col)) != 1
                  end
    return df[:,col_changed]
end

function delete_deleted_imgs()
    df = CSV.read("plots/df.csv", DataFrame)
    imgs = readdir("plots/new_plots")
    imgs = map(img -> img[1:end-4], imgs)
    idxs = map(img -> img in imgs, df.Image_Name)
    CSV.write("plots/df.csv", df[idxs,:])
end

# df = CSV.read("plots/df.csv", DataFrame)
# df = filter_unchanged_cols(df) 
# df_quench  = @subset(df, :Space_Time .== "quench")
# df_rindler = @subset(df, :Space_Time .== "rindler")
# df_flat    = @subset(df, :Space_Time .== "flat");

# filter_unchanged_cols(df_rindler)
delete_deleted_imgs()