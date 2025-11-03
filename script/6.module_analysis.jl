using RCall
using CSV, DataFrames

###########################################
# helper funtion to run  the WGCNA from R #
###########################################

Rwgcna <- function(cutoff) {
  merge_matrix <- read.csv("precomputed/wgcna_merge_matrix.csv",row.names = 1)
  geneSD = apply(merge_matrix, 1, sd)
  merge_matrix = merge_matrix[ geneSD > summary(geneSD)[2], colSums(merge_matrix) > 0]

  meta = read.csv("precomputed/wgcna_30bin_anno.csv")
  celltype = meta$celltype
  nP = meta$proteinNum

  input_mat = t(scale(merge_matrix)) # gene in column for wgcna

  powers = c(c(1:10), seq(from = 12, to = 30, by = 2))
  sft = pickSoftThreshold(
    input_mat ,             # <= Input data
    powerVector = powers,
    verbose = 5
  )

  picked_power = sft$powerEstimate
  picked_power = 4

  adjacency <- adjacency(input_mat,
                         power = picked_power,
                         type = "signed" )
  TOMadj <- TOMsimilarity(adjacency)

  dissTOMadj <- 1- TOMadj
  hclustGeneTree <- hclust(as.dist(dissTOMadj), method = "ward.D2")

  minModuleSize = 50
  dynamicMods <- cutreeDynamic(dendro = hclustGeneTree,
                               distM = dissTOMadj,
                               cutHeight = 1,
                               deepSplit = 2, # large => more module
                               pamRespectsDendro = F,
                               minClusterSize = minModuleSize)

  table(dynamicMods)
  dynamicColors <- labels2colors(dynamicMods)
  table(dynamicColors)


  dynamic_MEDissThres <- cutoff

  merge_dynamic_MEDs <- mergeCloseModules(input_mat, dynamicColors, cutHeight = dynamic_MEDissThres, verbose = 3)

  dynamic_mergedColors <- merge_dynamic_MEDs$colors

  mergedMEs <- merge_dynamic_MEDs$newMEs

  pdf(paste0("cutTree",dynamic_MEDissThres,".pdf"))

  plotDendroAndColors(hclustGeneTree, cbind(dynamicColors, dynamic_mergedColors),
                      c("Dynamic Tree Cut", "Merged dynamic"),
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  dev.off()

  moduleColors <- dynamic_mergedColors
  MEs <- mergedMEs

  MEs <- moduleEigengenes(input_mat, moduleColors)$eigengenes
  MEs <- orderMEs(MEs)

  module_df <- data.frame(
    gene_id = colnames(input_mat),
    colors =  merge_dynamic_MEDs$colors
  )

  MEs1 <- MEs[,setdiff(colnames(MEs),c("MEgrey"))]

  cor_ <- cor(MEs1,1:nrow(MEs1))
  p_ <- corPvalueStudent(cor_, nSample = dim(input_mat)[2])


  MEs0 <- orderMEs(MEs1)
  return(list(value = MEs0, df = module_df))
}


wgcna(cutoff) = begin
    @rput cutoff
    R""" res = Rwgcna(cutoff);MEs0 = res$value;module_df = res$df"""
    @rget module_df
    @rget MEs0
    ( module_df ,MEs0 )
end


########################
# get the wgcna result #
########################
module_df, module_value = wgcna(0.3)

# order the module by activity value
function moduleOrder(vec)
    a = sortperm(vec)[end]
end

a = Matrix(module_value)
colorLabel = names(module_value)
cor_ = abs.(corspearman(a,1:30)) |> vec
# colorLabel =  colorLabel[idx]
modod = mapslices( moduleOrder, a, dims = 1) |> vec

Ts = colorLabel[sortperm(modod)]
Ts = replace.(Ts, "ME" => "")
DataFrame( :color => Ts, :time => modod[sortperm(modod)], :cor => cor_ )

module_value = module_value[:,sortperm(modod)]
hm = module_value

meta = CSV.read("precomputed/wgcna_30bin_anno.csv",DataFrame)
celltype = meta.celltype

###################################
# plot the protein activity curve #
# fig 6. C                        #
###################################

pltdf = module_value
pltdf.time = 1:30

pltdf = stack(pltdf, Not(:time), :time)
pltdf.value .= Float64.(pltdf.value)

pltdf.col .= "1"
idx = pltdf.variable .∈ Ref([ "MEyellow", "MEblue", "MEtan"])
pltdf.col[idx] .= "2"

pltdf.celltype .= "IPC-EN"
pltdf.celltype[ pltdf.time .∈ Ref(1:7)] .= "RG"
pltdf.celltype[ pltdf.time .∈ Ref(21:30)] .= "EN"

plt = data(pltdf) * mapping(:time, :value, color = :celltype, row = :variable) *
    (visual(Scatter))

plt2 = data(pltdf) * mapping(:time, :value, row = :variable) * smooth(degree = 2);

fig = draw(plt + plt2; axis = (; width = 120, height = 120))


##############################################
# output the module protein for visulization #
##############################################
#+begin_src julia :eval no

R"""
    M = cor(t(merge_matrix))
    N = rownames(merge_matrix)
"""
@rget M
@rget N
module_df = CSV.read( "precomputed/module_df.csv", DataFrame)

# pick module Gene
df = DataFrame(M, :auto)
rename!(df, names(df) .= N)
color_ = ["yellow","blue"]
idx = N .∈ Ref( module_df.gene_id[ module_df.colors .∈ Ref(color_) ])

net = df[idx,idx]

source = []
target = []
weight = []
cutoff = 0.83

for i in 2:size(net,1)-1
    for j in  i+1:size(net,2)
        w = net[i,j]
        if w > cutoff
            push!(source, names(net)[i] )
            push!(target,  names(net)[j] )
            push!(weight,  w )
        end
    end
end

subNET = DataFrame( :source => source, :target => target, :weight => weight)
nodeInfo = DataFrame( :Id => vcat(source, target) |> unique )

nodeInfo.Label .= ""
nodeInfo.Label[ idx1  ] .= nodeInfo.Id[idx1]
nodeInfo.Label[ idx2  ] .= nodeInfo.Id[idx2]
leftjoin!(nodeInfo, module_df, on = :Id => :gene_id)
leftjoin!(subNET, module_df, on = :source => :gene_id)

CSV.write( "network.edge.tsv", subNET, delim = "\t" )
CSV.write( "network.node.tsv", nodeInfo, delim = "\t" )
