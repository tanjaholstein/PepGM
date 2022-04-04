
#configuration
WorkflowPath = os.path.dirname(os.path.realpath(workflow.snakefile))
include: WorkflowPath +'/snakefiles/config.snakefile'
include : WorkflowPath + '/snakefiles/searchgui+pepshaker.snakefile'
include : WorkflowPath + '/snakefiles/buildDatabase.snakefile'
include : WorkflowPath + '/snakefiles/filterHostCrap.snakefile'
#include: WorkflowPath + '/snakefiles/pepGM.snakefile'



rule all:
     input:
          expand(ResultsDir + SampleName + '/Prior{Prior}/'+ReferenceDBName +'_PepGM_Results_a{alpha}_b{beta}_p{Prior}.png',alpha = AlphaRange, beta = BetaRange, Prior = prior),
          ResultsDir +SampleName +'/paramcheck.png',
#          DataDirectory + SampleName + '/' + firstTarget +'_proteinCount.png'


rule getTargets:
     input:
          ResultsDir + SampleName+'/{DBname}_Default_PSM_Report.txt',
          ResourcesDir + TaxidMapping + 'accessions_hashed.npy',
          ResourcesDir + TaxidMapping + 'taxids.txt'
     params:
          query = ResultsDir + SampleName+'/{DBname}_query_accessions.txt',
          samplename = SampleName,
          hostname = HostName,
          DBname = ReferenceDBName
     conda: 'envs/graphenv.yml'
     output:
          ResultsDir +SampleName+'/{DBname}_mapped_taxids.txt'
     shell: "python3 workflow/scripts/getTargets.py -rq {input[0]} -q {params.query} -d {input[1]} -t {input[2]} -r {output} "


rule CreateFactorGraph:
     input: ResultsDir +SampleName+'/{DBname}_Default_PSM_Report.txt',
            ResultsDir + SampleName+'/{DBname}_mapped_taxids.txt'
     output: ResultsDir +SampleName+'/{DBname}_PepGM_graph.graphml'
     conda: 'envs/graphenv.yml'
     params:
          samplename = SampleName,
          hostname = HostName,
          DBname = ReferenceDBName,
          targetTaxa = TargetTaxa,
          sourceDB = sourceDB,
          ResultsDir = ResultsDir
     log: ResultsDir +SampleName+'/{DBname}.log'
     shell: 'cp config/config.yaml {params.ResultsDir}/{params.samplename}/  && python3 workflow/scripts/CreatePepGMGraph.py --sourceDB {params.sourceDB} --targetTaxa {input[1]}  --PSM_Report {input[0]} --PeptideMapPath '+ResultsDir+'{params.samplename}/Mapped_Taxa_Proteins.json --out {output} &>> {log}'


rule RunPepGM:
     input: ResultsDir +SampleName+'/{DBname}_PepGM_graph.graphml'
     output: ResultsDir +SampleName+'/Prior{Prior}/{DBname}_PepGM_Results_a{alpha}_b{beta}_p{Prior}.csv'
     conda: 'envs/graphenv.yml'
     params:
          samplename = SampleName,
          hostname = HostName,
          DBname = ReferenceDBName,
          targetTaxa = TargetTaxa,
          firstTarget = firstTarget
     log : ResultsDir + SampleName+'/Prior{Prior}/log/{DBname}_a{alpha}_b{beta}_pepGM.log'
     shell: 'python3  workflow/scripts/PepGM.py --GraphMLPath {input} --prior {wildcards.Prior} --alpha {wildcards.alpha} --beta {wildcards.beta} --out {output} &>> {log}'


rule BarPlotResults:
     input: ResultsDir +SampleName+'/Prior{Prior}/{DBname}_PepGM_Results_a{alpha}_b{beta}_p{Prior}.csv'
     output:  ResultsDir +SampleName+'/Prior{Prior}/{DBname}_PepGM_Results_a{alpha}_b{beta}_p{Prior}.png',
     conda: 'envs/graphenv.yml'
     params: NumberofResults = TaxaInPlot
     log : ResultsDir + SampleName+'/Prior{Prior}/log/{DBname}_a{alpha}_b{beta}.log'
     shell: 'python3 workflow/scripts/BarPlotResults.py --ResultsFile {input} --NumberofResults {params.NumberofResults} --out {output} &>> {log}'


rule FindBestParameters:
     input: 
          expand(ResultsDir + SampleName + '/Prior{Prior}/'+ ReferenceDBName +'_PepGM_Results_a{alpha}_b{beta}_p{Prior}.png',alpha = AlphaRange, beta = BetaRange, Prior = prior)
     output: ResultsDir +SampleName+'/paramcheck.png'
     conda: 'envs/graphenv.yml'
     log : ResultsDir + SampleName+'/paramcheck.log'
     params: Results = ResultsDir+SampleName,
             host = ScientificHostName
     shell: 'python3 workflow/scripts/GridSearchAnalysis.py --host {params.host} --resultsfolder {params.Results} --out {output} &>> {log}'


#rules only needed for debugging/checking results
rule CountProteinsInJson:
     input: ResultsDir +SampleName+'/{taxonname}.json'
     output : ResultsDir + SampleName+'/{taxonname}_proteinCount.png'
     params : NumberofResults = TaxaInProteinCount
     conda: 'envs/graphenv.yml'
     shell: 'python3 workflow/scripts/CountProteins.py --ResultsFile {input} --NumberofResults {params.NumberofResults} --out {output}'
