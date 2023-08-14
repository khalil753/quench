using CSV, DataFrames, DataFramesMeta

function filter_unchanged_cols(df)
    col_changed = map(eachcol(df)) do col
                      length(Set(col)) != 1
                  end
    return df[:,col_changed]
end

function delete_deleted_imgs(path)
    df = CSV.read("$path/df.csv", DataFrame)
    imgs = readdir("$path")
    imgs = map(img -> img[1:end-4], imgs)
    idxs = map(img -> img in imgs, df.Image_Name)
    CSV.write("$path/df.csv", df[idxs,:])
end
function delete_deleted_rows(path)
    df = CSV.read("$path/df.csv", DataFrame)
    files = readdir("$path")
    imgs = filter(img_name -> img_name[end-3:end]==".pdf", files)
    for img in imgs
        if img[1:end-4] in df.Image_Name continue end
        rm("$path/$img"); rm("$path/$(img[1:end-4]).jld")
    end
end

function sort_alfabetically(path)
  df = CSV.read( "$path/df.csv", DataFrame)
  println(df)
  sort!(df, [:Image_Name])
  println(df)
  CSV.write("$path/df.csv", df)    
end

# df = CSV.read("plots/df.csv", DataFrame)
# df = filter_unchanged_cols(df) 
# df_quench  = @subset(df, :Space_Time .== "quench")
# df_rindler = @subset(df, :Space_Time .== "rindler")
# df_flat    = @subset(df, :Space_Time .== "flat");

# filter_unchanged_cols(df_rindler)
delete_deleted_imgs("plots/new_plots/quench_2")
# imgs = delete_deleted_rows("plots/new_plots/quench")
# sort_alfabetically("plots/new_plots/quench")
