%include('header',title='confirm')
<head>
<meta name='viewport' content='width=device-width, initial-scale=1, maximum-scale=1'>
</head>
<body>
%include('navbar')
%include('apps/alert')
<div class="container-fluid">
<form class="form-horizontal" action="/confirm" method="post" novalidate>
<input type="hidden" name="app" value="{{app}}">
%if defined('cid'):
<input type="hidden" name="cid" value="{{cid}}">
%end

<a href="https://bio.cst.temple.edu/~hey/program_files/FPG/FPG_Documentation.htm" class="help btn btn-info" target="_blank"><span class="glyphicon glyphicon-question-sign"></span></a>

<div class="col-xs-12" style="height:5px"></div>
<div class="col-xs-12 visible-xs" style="height:5px"></div>
<div class="form-group">
	<div class="col-xs-2">
		<button type="submit" class="btn btn-success"> <!-- pull-right -->
		Continue <em class="glyphicon glyphicon-forward"></em> </button>
	</div>
	<label for="desc" style="text-align:right" class="control-label col-sm-4 hidden-xs">
		<a href="#" data-toggle="tooltip" title="Separate labels by commas">Labels:</a></label>
	<div class="hidden-xs col-sm-6">
		<input type="text" id="desc" name="desc" class="form-control" style="width:100%"
			data-role="tagsinput" title="e.g. v2.5.1,bottleneck">
	</div>
</div>
<div class="col-xs-12" style="height:5px"></div>

<ul class="nav nav-pills" role="tablist">
	<li role="presentation">
		<a href="#BASIC" aria-controls="home" role="tab" data-toggle="tab">BASIC</a>
	</li>
	<li role="presentation">
		<a href="#MUTATIONS" aria-controls="home" role="tab" data-toggle="tab">MUTATIONS</a>
	</li>
	<li role="presentation">
		<a href="#SELECTION" aria-controls="home" role="tab" data-toggle="tab">SELECTION</a>
	</li>
	<li role="presentation">
		<a href="#POPULATION" aria-controls="home" role="tab" data-toggle="tab">POPULATION</a>
	</li>
	<li role="presentation">
		<a href="#COMPUTATION" aria-controls="home" role="tab" data-toggle="tab">COMPUTATION</a>
	</li>
</ul>

<div class="tab-content">

<div role="tabpanel" class="tab-pane fade in active" id="BASIC">
	<div class="form-group">
		<label for="g_popsize" class="control-label col-xs-6">
			* popsize &mdash; the # of haploid genome copies (gametes) in the population:
		</label>
		<div class="col-xs-12 col-sm-6">
			<input type="number" class="form-control" name="g_popsize" value="{{g_popsize}}"/>
		</div>
	</div>
	<div class="form-group">
		<label for="o_polytype" class="control-label col-xs-6">
			* polytype &mdash; types of polymorphisms to analysis:<br> 
              &nbsp;&nbsp;&nbsp; <font size="-1">(B=beneficial mutations, H=harmful mutations, N=neutral mutations)</font></label>
		<div class="col-xs-12 col-sm-6">
			<input type="text" class="form-control" name="o_polytype" value="{{o_polytype}}"/>
		</div>
	</div>
	<div class="form-group">
		<label for="c_chromesegs" class="control-label col-xs-6">
			* chromesegs &mdash; the # of segments per chromosome<br>
			  32 bits per segment:</label>
		<div class="col-xs-12 col-sm-6">
			<input type="number" class="form-control" name="c_chromesegs" value="{{c_chromesegs}}"/>
		</div>
	</div>
	<div class="form-group">
		<label for="x_chromenum" class="control-label col-xs-6">
			* chromenum &mdash; number of chromosomes:</label>
		<div class="col-xs-12 col-sm-6">
			<input type="number" class="form-control" name="x_chromenum" value="{{x_chromenum}}"/>
		</div>
	</div>
	<div class="form-group">
		<label for="v_pop_u_s" class="control-label col-xs-6">
			* pop_u_s &mdash; population selected mutation rate:</label>
		<div class="col-xs-12 col-sm-6">
			<input type="number" class="form-control" name="v_pop_u_s" value="{{v_pop_u_s}}"/>
		</div>
	</div>
	<div class="form-group">
		<label for="u_pop_u_n" class="control-label col-xs-6">
			* pop_u_n &mdash; population neutral mutation rate:</label>
		<div class="col-xs-12 col-sm-6">
			<input type="number" class="form-control" name="u_pop_u_n" value="{{u_pop_u_n}}"/>
		</div>
	</div>
	<div class="form-group">
		<label for="r_poprecrate" class="control-label col-xs-6">
			* poprecrate &mdash; population recombination rate:<br>
              &nbsp;&nbsp;&nbsp; <font size="-1">(gamsize &times; recombination rate per chromosome)</font></label>
		<div class="col-xs-12 col-sm-6">
			<input type="number" class="form-control" name="r_poprecrate" value="{{r_poprecrate}}"/>
		</div>
	</div>
</div>

<div role="tabpanel" class="tab-pane fade in inactive" id="MUTATIONS">
	<div class="form-group">
		<label for="f_del_prop" class="control-label col-xs-6">
			* del_prop &mdash; proportion of selected mutations that are deleterious:</label>
		<div class="col-xs-12 col-sm-6">
			<input type="number" class="form-control" name="f_del_prop" value="{{f_del_prop}}"/>
		</div>
	</div>
	<div class="form-group">
		<label for="h_dominance" class="control-label col-xs-6">
			* dominance &mdash; dominance of selected mutations:</label>
		<div class="col-xs-12 col-sm-6">
			<input type="number" class="form-control" name="h_dominance" value="{{h_dominance}}"/>
		</div>
	</div>
	<div class="form-group">
		<label for="y_epistasis" class="control-label col-xs-6">
			epistasis &mdash; degree of epistatis among selected mutations:</label>
		<div class="col-xs-12 col-sm-6">
			<input type="number" class="form-control" name="y_epistasis" value="{{y_epistasis}}"/>
		</div>
	</div>
	<div class="form-group">
		<label for="w_fit_schema" class="control-label col-xs-6">
			w_fit_schema:</label>
		<div class="col-xs-12 col-sm-6">
			<select class="form-control" name="w_fit_schema">
			%opts = {'': 'none', 'A': 'additive', 'M': 'multiplicative', 'E': 'epistatic'}
			%for key, value in opts.iteritems():
				%if defined('w_fit_schema'):
					%if key == w_fit_schema:
						<option selected value="{{key}}">{{value}}</option>
					%else:
						<option value="{{key}}">{{value}}</option>
					%end
				%else:
					<option value="{{key}}">{{value}}</option>
				%end
			%end
			</select>
		</div>
	</div>
</div>

<div role="tabpanel" class="tab-pane fade in inactive" id="SELECTION">
	<div class="form-group">
		<label for="s_popscoeff" class="control-label col-xs-6">
			* popscoeff &mdash; population selection coefficient:<br>
              &nbsp;&nbsp;&nbsp;<font size="-1">(gamsize &times; selection coefficient)</font>
        </label>
		<div class="col-xs-12 col-sm-6">
			<input type="number" class="form-control" name="s_popscoeff" value="{{s_popscoeff}}"/>
		</div>
	</div>
</div>

<div role="tabpanel" class="tab-pane fade in inactive" id="POPULATION">
	<div class="form-group">
		<label for="k_numpops" class="control-label col-xs-6">
			* numpops &mdash; number of subpopulations:</label>
		<div class="col-xs-12 col-sm-6">
			<input type="number" class="form-control" name="k_numpops" value="{{k_numpops}}"/>
		</div>
	</div>
	<div class="form-group">
		<label for="m_mrate" class="control-label col-xs-6">
			mrate &mdash; migration rate out of each population per generation:</label>
		<div class="col-xs-12 col-sm-6">
		%if defined('m_mrate'):
			<input type="number" class="form-control" name="m_mrate" value="{{m_mrate}}"/>
		%else:
			<input type="number" class="form-control" name="m_mrate"/>
		%end
		</div>
	</div>
</div>

<div role="tabpanel" class="tab-pane fade in inactive" id="COMPUTATION">
	<div class="form-group">
		<label for="n_measurements" class="control-label col-xs-6">
			* measurements &mdash; the number of measurements to make:</label>
		<div class="col-xs-12 col-sm-6">
			<input type="number" class="form-control" name="n_measurements" value="{{n_measurements}}"/>
		</div>
	</div>
	<div class="form-group">
		<label for="i_t_start" class="control-label col-xs-6">
			* t_start &mdash; length of burn-in period in units of population size generations:</label>
		<div class="col-xs-12 col-sm-6">
			<input type="number" class="form-control" name="i_t_start" value="{{i_t_start}}"/>
		</div>
	</div>
	<div class="form-group">
		<label for="j_t_between" class="control-label col-xs-6">
			t_between &mdash; time between measurements after burnin<br>
			in units of population size generations:</label>
		<div class="col-xs-12 col-sm-6">
			<input type="number" class="form-control" name="j_t_between" 
                        %if defined('j_t_between'):
                          value="{{j_t_between}}"
                        %end
                        />
		</div>
	</div>
	<div class="form-group">
		<label for="l_ld_analysis" class="control-label col-xs-6">
			LD_analysis &mdash; invokes linkage disequilibrium analysis: <br>
			<font size="-1">(some combination of R, S, D, P or B &mdash; see manual)</font></label>
		<div class="col-xs-12 col-sm-6">
		%if defined('l_ld_analysis'):
			<input type="text" class="form-control" name="l_ld_analysis" value="{{l_ld_analysis}}"/>
		%else:
			<input type="text" class="form-control" name="l_ld_analysis"/>
		%end
		</div>
	</div>
	<div class="form-group">
		<label for="t_pseudo_data" class="control-label col-xs-6">
			pseudo_data &mdash; causes a pseudo data set, in SITES Program format to be included in output:</label>
		<div class="col-xs-12 col-sm-6">
			<input type="checkbox" name="t_pseudo_data" value="true"
				%if t_pseudo_data== 'true':
				checked
				%end
			/>
		</div>
	</div>
	<div class="form-group">
		<label for="q_absolutefit" class="control-label col-xs-6">
			absolutefit &mdash; causes fitness measurements to include all those mutations that fixed:</label>
		<div class="col-xs-12 col-sm-6">
			<input type="checkbox" name="q_absolutefit" value="true"
				%if q_absolutefit== 'true':
				checked
				%end
			/>
		</div>
	</div>
</div>

</div>

* required parameter

</form>
%include('footer')
