# all
all: sn.sim10.stk

# don't delete intermediate ("secondary") targets
.SECONDARY:

# DART logfile-related functions
# the only one used in the rest of the file is RunAndLog
# format is RunAndLog(ProgramName,OutputStem,OutputSuffix,Args)
# this will create files '$*.OutputStem.OutputSuffix' (stdout) and '$*.OutputStem.log' (logfile)
# it also pipes the logfile through dartlog.pl, with nice colorful options
LogDirectives = -log 5
DartLogViewer = dartlog.pl -fancy -cat -min 6 -show ECFGDP
RemoveLog = (test -e $*.$(1).log && rm $*.$(1).log) || true
RunAndLog = @$(call RemoveLog,$(2)); \
            echo "$*.$(2).log : $(1) $(4) >$*.$(2).$(3)"; \
            tcsh -c '($(1) $(LogDirectives) -logtime -logcopy $*.$(2).log $(4) >$*.$(2).$(3)) |& $(DartLogViewer)'


# special rule for ../sn.eg and other codon grammars
../%.eg: %.eg
	xrate -g $< -x $*.tmp.eg dummy.stock -noa >/dev/null
	cat $*.tmp.eg | perl -ne 's/HYBRID/COD/g;s/_ancestral//g;s/_DEFAULT//g;if(/end chain/){$$n=1}elsif(/end hybrid/){$$n=0}print unless/_(root|x|y)/ || $$n'  >$@
	rm $*.tmp.eg

# build grammar from generator script & template
%.eg: %.tmpl.eg %.pl
	perl $*.pl $*.tmpl.eg >$*.eg

# expand grammar macros from Stockholm file
%.expanded.eg: %.eg %.stk
	$(call RunAndLog,xrate,expanded,out,$*.stk -g $*.eg -x $*.expanded.eg -noa)

# train grammar on Stockholm file
%.trained.eg: %.eg %.stk
	$(call RunAndLog,xrate,trained,out,$*.stk -g $*.eg -t $*.trained.eg -noa)

# extract tree(s) from Stockholm file
%.nh: %.stk
	cat $< | perl -ne 'if (/#=GF\s+NH\s+(\S.*)/) { print "$$1\n" }' >$@

# simulate from extracted tree & trained grammar
# 1 alignment
%.sim1.stk: %.trained.eg %.nh
	$(call RunAndLog,simgram,sim1,stk,-rndtime -g $*.trained.eg -t $*.nh -n 1)

# 10 alignments
%.sim10.stk: %.trained.eg %.nh
	$(call RunAndLog,simgram,sim10,stk,-rndtime -g $*.trained.eg -t $*.nh -n 10)

# 100 alignments
%.sim100.stk: %.trained.eg %.nh
	$(call RunAndLog,simgram,sim100,stk,-rndtime -g $*.trained.eg $*.nh -n 100)

# re-impute parameters from simulated data
# 1 alignment
%.imp1.eg: %.eg %.sim1.stk
	$(call RunAndLog,xrate,imp1,out,$*.sim1.stk -g $*.eg -t $*.imp1.eg -noa)

# 10 alignments
%.imp10.eg: %.eg %.sim10.stk
	$(call RunAndLog,xrate,imp10,out,$*.sim10.stk -g $*.eg -t $*.imp10.eg -noa)

# 100 alignments
%.imp100.eg: %.eg %.sim100.stk
	$(call RunAndLog,xrate,imp100,out,$*.sim100.stk -g $*.eg -t $*.imp100.eg -noa)

# compare trained to simulated-imputed parameters (1- & 10-alignment simulated datasets)
%.cmp.sxpr: %.trained.eg %.imp1.eg %.imp10.eg
	compareParams.pl $*.trained.eg $*.imp1.eg $*.imp10.eg >$@

# compare trained to simulated-imputed parameters (1-alignment simulated dataset only)
%.cmp1.sxpr: %.trained.eg %.imp1.eg
	compareParams.pl $*.trained.eg $*.imp1.eg >$@

# compare self to noself
%.noself.cmp.sxpr: %.trained.eg noself-%.trained.eg
	compareParams.pl $*.trained.eg noself-$*.trained.eg >$@

# digest of parameter comparison
# extracts (param ...) , ((def ...)) & ((obs ...)) clauses
# prints digest as columns:
# PARAM DEF1 DEF2 [DEF3 ....] OBS1 OBS2 [OBS3 ...]
%.digest: %.sxpr
	perl -e 'use DartSexpr;DartSexpr->from_file(shift)->map_child(sub{$$c=shift;print join(" ",$$c->param->value,map($$_->def->value,$$c->values),map($$_->obs->value,$$c->values)),"\n"})' $< >$@

# xgraph of parameter value comparison
%.xgdef: %.cmp.digest
	cat $< | fields 1 3 | xgraph -nl -m

# xgraph of parameter counts comparison
%.xgobs: %.cmp.digest
	cat $< | fields 4 6 | xgraph -nl -m


# 100 repetitions of %.cmp1.digest in separate subdirectories
%.100reps: %.eg %.trained.eg %.nh
	test -e $@ || mkdir $@
	perl -e '$$cwd=`pwd`;chomp$$cwd;for$$n(1..100){$$d="$@/$$n";mkdir$$d;for$$f(qw($*.eg $*.trained.eg $*.nh)){symlink"$$cwd/$$f","$$d/$$f"};system"make $$d/$*.cmp1.digest"}'

# scatterplot of k-values from 100 reps
# example of analysis:
# cat sn.k-xyplot | perl -e 'while(<>){($x,$y)=split;++$n{$x};$t{$x}+=$y;$tt{$x}+=$y*$y}for$x(sort keys%n){$m=$t{$x}/$n{$x};$v=$tt{$x}/$n{$x}-$m*$m;print"$x $m ",sqrt($v),"\n"}'
%.k-xyplot: %.100reps
	egrep "^k_" $*.100reps/*/$*.cmp1.digest | fields 1 2 >$@

# copy data from sn model to other models
hky85.stk noself-sn.stk: sn.stk
	ln -sf $< $@

# dummy, empty files that need to be there to trigger make rules...
# various grammar templates
ka_ks.tmpl.eg phylohmmgenerator.tmpl.eg:
	touch $@

# "noself" versions of grammar templates
noself-%.tmpl.eg:
	touch $@

# "cat" perl scripts
hky85.pl: cat.pl
	ln -sf $< $@

# "noself" versions of perl scripts
noself-%.pl: noself.pl
	ln -sf $< $@
