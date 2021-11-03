""" Benchmarking the OptiFit algorithm using datasets split into reference & query sets """

configfile: 'config/config.yaml'

printrefs = ['t'] # include the reference sequences in the final MCC score for the most fair comparison to OptiClust.
dist_thresh = config['dist_thresh']
datasets = config['datasets']
methods = config['methods']
seeds = range(1, config['seeds'] + 1)
random_refweight_options = config['weights']

fracs = config['ref_fracs']
ref_fracs = [i/10 for i in range(fracs['start'], fracs['stop'], fracs['step'])]

weights_to_fracs = {'simple': ref_fracs,
                    'abundance': [0.5],
                    'distance': [0.5]
}
for weight in set(weights_to_fracs.keys()) - set(random_refweight_options):
    weights_to_fracs.pop(weight)

results_files = [f'results/{dataset}/refweight_{ref_weight}/reffrac_{ref_frac}/samplefrac_{round(1 - ref_frac, 1)}/seed_{seed}/fit/method_{method}/printref_{printref}/results.tsv'
    for dataset in datasets
    for ref_weight, ref_fracs in weights_to_fracs.items()
    for ref_frac in ref_fracs
    for seed in seeds
    for method in methods
    for printref in printrefs]

onstart:
    print(f"🎛  Params will produce {len(results_files)} OptiFit runs")

wildcard_constraints:
    seed="\d+",
    printref="t|f"

rule rbind_optifit_split:
    input:
        code='code/py/rbind_files.py',
        tsv=results_files
    output:
        tsv='results/optifit_split_results.tsv'
    log:
        'log/rbind_optifit_split.txt'
    script:
        'code/py/rbind_files.py'

rule render_readme:
    input:
        rmd='README.Rmd'
    output:
        md="README.md"
    shell:
        """
        R -e "rmarkdown::render('{input.rmd}')"
        """

rule split_weighted_subsample:
    input:
        code="code/py/split_weighted_subsample.py",
        fasta="data/{dataset}/processed/{dataset}.fasta",
        count="data/{dataset}/processed/{dataset}.count_table",
        dist="data/{dataset}/processed/{dataset}.dist"
    output:
        ref_accnos="data/{dataset}/refweight_{ref_weight}/reffrac_{ref_frac}/samplefrac_{sample_frac}/seed_{seed}/ref.accnos",
        query_accnos="data/{dataset}/refweight_{ref_weight}/reffrac_{ref_frac}/samplefrac_{sample_frac}/seed_{seed}/sample.accnos",
        all_accnos="data/{dataset}/refweight_{ref_weight}/reffrac_{ref_frac}/samplefrac_{sample_frac}/seed_{seed}/combined/{dataset}.accnos"
    log:
        "log/{dataset}/split_weighted_subsample.refweight_{ref_weight}.reffrac_{ref_frac}.samplefrac_{sample_frac}.seed_{seed}.txt"
    params:
        dissim_thresh=dist_thresh
    script:
        "code/py/split_weighted_subsample.py"

rule check_split:
    input:
        code="code/py/check_split.py",
        ref=rules.split_weighted_subsample.output.ref_accnos,
        query=rules.split_weighted_subsample.output.query_accnos,
        all=rules.split_weighted_subsample.output.all_accnos
    output:
        txt=temp("data/{dataset}/refweight_{ref_weight}/reffrac_{ref_frac}/samplefrac_{sample_frac}/seed_{seed}/check_split.txt")
    script:
        'code/py/check_split.py'

rule prep_subsample_ref:
    input:
        accnos=rules.split_weighted_subsample.output.ref_accnos,
        count="data/{dataset}/processed/{dataset}.count_table",
        fasta="data/{dataset}/processed/{dataset}.fasta",
        dist="data/{dataset}/processed/{dataset}.dist",
        mothur="bin/mothur-1.47.0/mothur"
    output:
        fasta=temp("data/{dataset}/refweight_{ref_weight}/reffrac_{ref_frac}/samplefrac_{sample_frac}/seed_{seed}/ref/{dataset}.pick.fasta"),
        count=temp("data/{dataset}/refweight_{ref_weight}/reffrac_{ref_frac}/samplefrac_{sample_frac}/seed_{seed}/ref/{dataset}.pick.count_table"),
        dist=temp("data/{dataset}/refweight_{ref_weight}/reffrac_{ref_frac}/samplefrac_{sample_frac}/seed_{seed}/ref/{dataset}.pick.dist")
    log:
        "log/{dataset}/prep_subsample_ref.refweight_{ref_weight}.reffrac_{ref_frac}.samplefrac_{sample_frac}.seed_{seed}.txt"
    params:
        outdir="data/{dataset}/refweight_{ref_weight}/reffrac_{ref_frac}/samplefrac_{sample_frac}/seed_{seed}/ref/"
    shell:
        """
        {input.mothur} "#set.logfile(name={log}); 
            set.dir(output={params.outdir});
            get.seqs(accnos={input.accnos}, fasta={input.fasta});
            get.seqs(accnos=current, count={input.count});
            get.dists(accnos=current, column={input.dist}) 
            "
        """

rule cluster_ref:
    input:
        dist=rules.prep_subsample_ref.output.dist,
        count=rules.prep_subsample_ref.output.count,
        mothur="bin/mothur-1.47.0/mothur"
    output:
        list=temp("results/{dataset}/refweight_{ref_weight}/reffrac_{ref_frac}/samplefrac_{sample_frac}/seed_{seed}/cluster/{dataset}.pick.opti_mcc.list"),
        sensspec=temp('results/{dataset}/refweight_{ref_weight}/reffrac_{ref_frac}/samplefrac_{sample_frac}/seed_{seed}/cluster/{dataset}.pick.opti_mcc.sensspec')
    params:
        outdir="results/{dataset}/refweight_{ref_weight}/reffrac_{ref_frac}/samplefrac_{sample_frac}/seed_{seed}/cluster/"
    log:
        "log/{dataset}/cluster.refweight_{ref_weight}.reffrac_{ref_frac}.samplefrac_{sample_frac}.seed_{seed}.log"
    benchmark:
        "benchmarks/{dataset}/cluster.refweight_{ref_weight}.reffrac_{ref_frac}.samplefrac_{sample_frac}.seed_{seed}.txt"
    resources:
        procs=8
    priority: 1
    shell:
        """
        {input.mothur} "#set.logfile(name={log}); 
            set.dir(output={params.outdir});
            set.seed(seed={wildcards.seed});
            set.current(processors={resources.procs});
            cluster(column={input.dist}, count={input.count}, cutoff=0.03) 
            "
        """

rule fit_subsample:
    input:
        fasta="data/{dataset}/processed/{dataset}.fasta",
        count="data/{dataset}/processed/{dataset}.count_table",
        column="data/{dataset}/processed/{dataset}.dist",
        reflist=rules.cluster_ref.output.list,
        mothur="bin/mothur-1.47.0/mothur"
    output:
        sensspec=temp('results/{dataset}/refweight_{ref_weight}/reffrac_{ref_frac}/samplefrac_{sample_frac}/seed_{seed}/fit/method_{method}/printref_{printref}/{dataset}.optifit_mcc.sensspec'),
        list=temp('results/{dataset}/refweight_{ref_weight}/reffrac_{ref_frac}/samplefrac_{sample_frac}/seed_{seed}/fit/method_{method}/printref_{printref}/{dataset}.optifit_mcc.list'),
        sum=temp('results/{dataset}/refweight_{ref_weight}/reffrac_{ref_frac}/samplefrac_{sample_frac}/seed_{seed}/fit/method_{method}/printref_{printref}/current_files.summary')
    params:
        outdir="results/{dataset}/refweight_{ref_weight}/reffrac_{ref_frac}/samplefrac_{sample_frac}/seed_{seed}/fit/method_{method}/printref_{printref}/"
    log:
        "log/{dataset}/fit.method_{method}.printref_{printref}.refweight_{ref_weight}.reffrac_{ref_frac}.samplefrac_{sample_frac}.seed_{seed}.txt"
    benchmark:
        "benchmarks/{dataset}/fit.method_{method}.printref_{printref}.refweight_{ref_weight}.reffrac_{ref_frac}.samplefrac_{sample_frac}.seed_{seed}.txt"
    resources:
        procs=8
    priority: 2
    shell:
        """
        {input.mothur} "#set.logfile(name={log}); 
            set.dir(output={params.outdir});
            set.seed(seed={wildcards.seed});
            set.current(processors={resources.procs});
            cluster.fit(reflist={input.reflist}, fasta={input.fasta}, count={input.count}, column={input.column},  printref={wildcards.printref}, method={wildcards.method}) 
            "
        """

rule list_seqs:
    input:
        list=rules.fit_subsample.output.list,
        mothur="bin/mothur-1.47.0/mothur"
    output:
        list_accnos=temp('results/{dataset}/refweight_{ref_weight}/reffrac_{ref_frac}/samplefrac_{sample_frac}/seed_{seed}/fit/method_{method}/printref_{printref}/{dataset}.optifit_mcc.accnos')
    params:
        outdir="results/{dataset}/refweight_{ref_weight}/reffrac_{ref_frac}/samplefrac_{sample_frac}/seed_{seed}/fit/method_{method}/printref_{printref}/"
    log:
        "log/{dataset}/list_seqs.method_{method}.printref_{printref}.refweight_{ref_weight}.reffrac_{ref_frac}.samplefrac_{sample_frac}.seed_{seed}.txt"
    priority: 3
    shell:
        """
        {input.mothur} "#set.logfile(name={log}); 
            set.dir(output={params.outdir});
            list.seqs(list={input.list})
            "
        """

rule count_input_sizes:
    input:
        py='code/py/count_input_sizes.py',
        fasta="data/{dataset}/processed/{dataset}.fasta",
        ref_accnos=rules.split_weighted_subsample.output.ref_accnos,
        query_accnos=rules.split_weighted_subsample.output.query_accnos
    output:
        tsv=temp("results/{dataset}/refweight_{ref_weight}/reffrac_{ref_frac}/samplefrac_{sample_frac}/seed_{seed}/input_size.tsv")
    script:
        'code/py/count_input_sizes.py'

rule fraction_reads_mapped3:
    input:
        code="code/py/fraction_reads_mapped.py",
        mapped=rules.list_seqs.output.list_accnos,
        query=rules.split_weighted_subsample.output.query_accnos,
        ref=rules.split_weighted_subsample.output.ref_accnos
    output:
        txt=temp('results/{dataset}/refweight_{ref_weight}/reffrac_{ref_frac}/samplefrac_{sample_frac}/seed_{seed}/fit/method_{method}/printref_{printref}/fraction_reads_mapped.txt')
    priority: 4
    script:
        "code/py/fraction_reads_mapped.py"

rule cbind_optifit_seed3:
    input:
        tsv=[rules.fit_subsample.output.sensspec,
            rules.fit_subsample.benchmark,
            rules.count_input_sizes.output.tsv,
            rules.fraction_reads_mapped3.output.txt,
            rules.check_split.output.txt],
        code='code/R/cbind_optifit_seed.R',
        fcns='code/R/functions.R'
    output:
        tsv='results/{dataset}/refweight_{ref_weight}/reffrac_{ref_frac}/samplefrac_{sample_frac}/seed_{seed}/fit/method_{method}/printref_{printref}/results.tsv'
    priority: 5
    script:
        'code/R/cbind_optifit_seed.R'

