using CSV, DataFrames, DataFramesMeta, Statistics, Plots, StatsPlots

AllTissues = CSV.read("AllTissues_Accessions.csv", DataFrame)

@linq AllTissues |> select!(:Accession, :Description, :ShootTip, :RootTip, :Callus, :Bud, :Xylem, :Leaf, :Bark)


function FindGene(df, name)
	findings = []
	for row in eachrow(df)
		if occursin(name,row.Description)
			push!(findings, row.Description)
		end
	end
	return(findings)
end


ShootSpecific = DataFrame(Accession = AllTissues.Accession, 
						  Description = AllTissues.Description,
						  RootTip = AllTissues.ShootTip./(AllTissues.RootTip.+1),
						  Bud = AllTissues.ShootTip./(AllTissues.Bud.+1),
						  Leaf = AllTissues.ShootTip./(AllTissues.Leaf.+1),
						  Xylem = AllTissues.ShootTip./(AllTissues.Xylem.+1),
						  Bark = AllTissues.ShootTip./(AllTissues.Bark.+1),
						  Callus = AllTissues.ShootTip./(AllTissues.Callus.+1))

ShootSpecific.Sum = sum.(eachrow(ShootSpecific[:,3:end]))
sort!(ShootSpecific, :Sum, rev = true)
  
@df ShootSpecific plot(1:nrow(ShootSpecific),:Sum, grid = false, label = "ShootTip", xlabel = "Genes", ylabel = "ShootTip frequency Over All Tissue")

Wuschels = FindGene(ShootSpecific, "WUSCHEL")

Wuschels_lacations = Dict()
for wus in Wuschels
    push!(Wuschels_lacations, "$wus" => findfirst(ShootSpecific.Description .== wus))
end


PLTs = FindGene(ShootSpecific,"PLT")
ESRs = FindGene(ShootSpecific,"ESR")
# WIND1 not found


Gene_list = [Wuschels; ESRs; PLTs]

Gene_lacations = DataFrame(Name = [], Location = [])
for gene in Gene_list
    push!(Gene_lacations, [gene findfirst(ShootSpecific.Description .== FindGene(ShootSpecific,gene)[1])])
end

sort!(Gene_lacations, :Location)
CSV.write("Gene_lacations.csv", Gene_lacations)