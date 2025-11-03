# including("ORBIT.jl") to load below data and function

pp, #theoretical peptide of protein
cds, #cds length
aa_len, #amino acid length
gravy, # GRAVY_Score
TM, # TM: TMhelice
ORBIT = let
	# theoretical peptide
	df = CSV.read("data/theoretical_peptide_counts.csv",DataFrame)
	pp = combine(groupby(df,:Gene), "Theoretical Peptide Count" => maximum => "n_peptide")
	# cds length
	cds = CSV.read("data/cds.txt",DataFrame)
	# protein length
	aa_len =  CSV.read("data/calculate_aminoacid_result.csv",DataFrame)
	# gravy
	df = CSV.read("data/gravy_results.csv",DataFrame)
	gravy = combine(groupby(df,:Gene), "GRAVY_Score" => maximum => "GRAVY_Score")
	# TMhmm
	df = CSV.read("data/TMhmm.txt",DataFrame)
	df.n = df.var"Transmembrane helices end" .= df.var"Transmembrane helices start"
	rename!(df, "Gene name" => "Gene")
	df = combine(groupby(df,:Gene), "n" => maximum => "n_TMhelice")
	df.n_TMhelice[ ismissing.(df.n_TMhelice)] .= 0
	TM = dropmissing(df)

	ORBIT(df::DataFrame;z_score::Bool=true,r = 0.5) = begin
		pro = names(df)[1]
		cols = names(df)[2:end]
		# add biochemistry info
		df = innerjoin(df, pp, on = Symbol(pro) => :Gene)
		df = innerjoin(df, gravy, on = Symbol(pro) => :Gene)
		df = innerjoin(df, TM, on = Symbol(pro) => :Gene)
		# new df for output
		newdf = similar(df[!,1:end-3])
		newdf[:,2:end] .= 0
		rename!(newdf, names(newdf) .=> names(df)[1:end-3])
		newdf[:,1] .= df[:,1]
		# regression adjustment
		for c in cols
			# should only correct those non-zero data.
			idx = df[:,c] .> 0
			X = Matrix(df[idx,end-2:end])
			y = Matrix(df[idx,[c]])

			F = ridge(X, y, r,bias=true)
			# extract regression weight
			A, b = F[1:end-1,:], F[end,:]
			newdf[idx,[c]] .= y .- X * A
		end
		# sort by the protein name
		newdf = sort(newdf,pro)
		# zscore
		if z_score
			newdf[:,2:end] = mapslices(zscore,Matrix(newdf[:,2:end]),dims=1)
		end
		return newdf
	end

	# export clean data and function
	pp, cds, aa_len, gravy, TM, ORBIT
end
