clear
set more off
set scheme uncluttered
scalar drop _all
capture log close
log using CDLT.log, replace
set type double

* --------------------------------
* TITLE: CDLT.do
* CREATED BY: Elizabeth Wrigley-Field
* CREATED: June 3, 2020
* LAST CHANGED: July 1, 2020
* --------------------------------

/* 
Contents: Make cause-deleted life table of police violence based on
data and code from Edwards, Lee, and Esposito (ELE) 2019.

Datasets needed:
	* fe_post_tables.dta (ELE)
	* fe_imp.dta (ELE)
	* Excel files Table_* : US 2017 lifes 
		(from: https://www.cdc.gov/nchs/data/nvsr/nvsr68/nvsr68_07-508.pdf )
	* popsize_2016_2011.dta : US population data from ACS, via IPUMS-USA

The ELE replication package is here:
https://osf.io/c8qxh/ 

The data files called up here that begin fe_* were made using these 
ELE files in R:
	model_ver_fin
	lifetables
	lifetable_multiple
	figures_tables

From the ELE replication package downloaded June 2, 2020, I 
made these changes: 
- model_ver_FIN.r, line 26: changed fe_pop_imputed_08_18 to fe_pop_imputed_13_18
- lifetable.R, line 8: changed ungroup(nat_dat) to ungroup
- lifetable.R, lines 43-44: left the file as-is but altered the estimation
of variables "L" (life table nLx) and "a" (life table nax) in the file below
- exported resulting datasets to Stata *.dta format
*/

* --------------------------------
* 0. Create nax values and true life table rates
* --------------------------------

local LTpath "NCHS-lifetables-2017"

/* 2017 Life Tables from the NCHS:
Table 02: All males
Table 03: All females
Table 11: Hispanic males
Table 12: Hispanic females
Table 14: non-Hispanic white males
Table 15: non-Hispanic white females
Table 17: non-Hispanic black males
Table 18: non-Hispanic black females
*/

local r02 Total
local r03 Total
local r11 Hispanic
local r12 Hispanic
local r14 White
local r15 White
local r17 Black
local r18 Black
local g02 Male
local g03 Female
local g11 Male
local g12 Female
local g14 Male
local g15 Female
local g17 Male
local g18 Female

foreach n in 02 03 11 12 14 15 17 18 {
		import excel "`LTpath'/Table`n'.xlsx",  cellrange(A3:G104) firstrow
		gen race = "`r`n''"
		gen gender = "`g`n''"
		save "`LTpath'/Table`n'.dta", replace
		clear
}
foreach n in 02 03 11 12 14 15 17 18 {
	append using "`LTpath'/Table`n'.dta"
}
tab race gender
ren gender sex // matches ELE variable name

* clean age variable
replace A = "100" if A=="100 and over"
split A, parse("–") destring
ren A1 age
tab age
drop A A2

* categorical race and gender variables
ren race racename
gen race = 1 if racename=="Black"
	replace race = 4 if racename=="Hispanic"
	replace race = 5 if racename=="White"
	replace race = 6 if racename=="Total"
	label define racelbl 1 "Black" 2 "Native" 3 "Asian/PI" ///
		4 "Latinx" 5 "White" 6 "Total"
	label values race racelbl

* Generate nax values for single-year observations
sort race sex age
gen nax = (Lx - lx[_n+1])/dx if age<100 
	// single-year obs, so no need to multiple lx by n (n=1)
replace nax = Lx/dx if age==100 // everyone dies in this category
sum nax, d
sum nax if age>0, d // .5 for all is very reasonable approximation

save "`LTpath'/lifetables.dta", replace


* Create abridged life table with 5-year exposures
clonevar ageorig = age
replace age = 1 if age<5 & age>1
replace age = 5*floor(age/5) if age>5 & age<85
replace age = 85 if age>85

collapse (sum) Lx_lifetable=Lx dx_lifetable=dx (max) lx_lifetable=lx, ///
	by(race sex age)
	// lx is max(lx) because that gives the lx for the beginning of each interval
gen mx_lifetable = dx/Lx
save Lxdx.dta, replace

** Make nax file ready to merge with ELE files
* Create true average age at death in single-year observations
* True average age at death is age at the bottom of the interval plus nax
use "`LTpath'/lifetables.dta", clear
gen agedeath = age + nax

* Collapse into 5-year abridged life table with open interval at 85
* For ages below 85, true age at death is dx-weighted average trueage
* within each interval
* For age 85, use e(85) as calculated by the NCHS and verify it equals
* nax as calculated by me
clonevar ageorig = age
replace age = 1 if age<5 & age>1
replace age = 5*floor(age/5) if age>5 & age<85
replace age = 85 if age>85

gen e85 = ex if ageorig==85

collapse (mean) agedeath=agedeath e85 [iweight=dx], by(race sex age)
gen nax = agedeath - age
table age, c(min nax max nax p50 nax)
table race sex if age==85, c(p50 nax p50 e85)
	* These should be equal, and are to at least 5 decimal places. --EWF, 6/7/2020
drop e85

* Create categories for Asian/Pacific Islander and Native
* Borrow API from white (API life expectancy is highest, white second-highets)
* and borrow Native from Black through age 50, Latinx afterward 
* (chosen because they are the closest mx matches)
egen tempgroup = group(sex age)
expandby 2, by(tempgroup) gen(fake1) // add Native observation
replace nax = . if fake1==1
replace race = 2 if fake1==1
sort sex age race
by sex age: egen tempnax1a = max(nax) if fake1==1 | race==1
by sex age: egen tempnax1b = max(nax) if fake1==1 | race==4
replace nax = tempnax1a if fake1==1 & age<=50
replace nax = tempnax1b if fake1==1 & age>50

expandby 2, by(tempgroup) gen(fake2) // add API observation
replace nax = . if fake2==1
replace race = 3 if fake2==1
sort sex age race
by sex age: egen tempnax2 = max(nax) if fake2==1 | race==5
replace nax = tempnax2 if fake2==1

drop temp* fake*
table race sex, c(p50 nax mean nax min nax max nax)

save nax.dta, replace

*** Estimate growth rate of population aged 85+ in specific race/sex groups,
* to estimate e(85) for Native and API groups for whom life tables aren't
* available
use popsize_2016_2011.dta, clear
	* Population sizes based on 5-year ACS in 2018 (centered 2016) and
	* 2013 (centered 2011)
reshape wide pop, i(race sex age) j(year)
gen growth = pop2018 / pop2013
gen r = ln(growth) / 5 // 5-year difference between 2011 and 2016
ren sex gender
gen sex = "Female" if gender==2
replace sex = "Male" if gender==1
foreach sex in Male Female {
	sum r if age==85 & sex=="`sex'" & race==2 // Native
	scalar r85_Nat_`sex' = r(mean)
	sum r if age==85 & sex=="`sex'" & race==3 // A & PI
	scalar r85_API_`sex' = r(mean)
}

* --------------------------------
* 1. Cause-deleted life tables for excessive force only
* --------------------------------
use "fe_post_tables.dta", clear
* use fe_median.dta, clear
	/*
	* only if this file --> alternative to fe_post_tables.dta
	ren race racename
	gen race = 1 if racename=="African American"
	replace race = 2 if racename=="American Indian/AK Native"
	replace race = 3 if racename=="Asian/Pacific Islander"
	replace race = 4 if racename=="Latinx"
	replace race = 5 if racename=="White"
	replace race = 6 if racename=="Total"
	*/

/* Notes on ELE notation
	q_med and m_med are the model-based estimates of deaths to policing,
		and m and dx are total mortality and total deaths
	L is nLx
*/	

* Merge in life table Lx values
merge 1:1 age sex race using Lxdx.dta
drop _merge
	
* Merge in nax values
merge m:1 age sex race using nax.dta
drop _merge
sort sex race age	

ren d_x dx_ELE
ren dx_lifetable dx
replace dx = dx_ELE if missing(dx)

clonevar lx_ELE = lx
replace lx = lx_lifetable if !missing(lx_lifetable)

* Correct the ELE L estimate, which is based on a model nax
* (implies mortality of .19 for all race/sex groups)
* My correction uses nax derived from NCHS life table for available
* race groups: Black, white, Latinx, and Total
* For Native and API, I use the Horiuchi and Coale (1982) formula (see PHG)
ren L Lx
replace Lx = Lx_lifetable if !missing(Lx_lifetable)
replace Lx = lx * nax if age==85 & inlist(race,1,4,5,6) // NCHS estimate

clonevar m_ELE = m
replace m = mx_lifetable if !missing(mx_lifetable)

local r2 Nat
local r3 API
foreach race in 2 3 {
	foreach sex in Male Female {
		replace Lx = ///
			lx * (1/m) * exp(-.0951 * r85_`r`race''_`sex' * m^(-1.4)) if ///
			age==85 & race==`race' & sex=="`sex'"
				// Horiuchi and Coale formula
	}
}

* Corrected life expectancy estimates using corrected Lx
ren e_med e_ELE
gen Tx_corrected = Lx if age==85
sort sex race age
forvalues a=80(-5)5 {
	replace Tx_corrected = Lx + Tx_corrected[_n+1] if age==`a'
}
foreach a in 1 0 {
	replace Tx_corrected = Lx + Tx_corrected[_n+1] if age==`a'
}
gen ex_corrected = Tx/lx
list age dx q_med m m_med lx Lx Tx ex if race==1 & sex=="Male"

* Compare my life expectancies to ELE's. Note nax fictional for Native & API
table race sex if age==85, c(p50 nax p50 ex_corrected p50 e_ELE)

* Check lifetime deaths to officer force; recreates main result of ELE paper
gen dxi = lx * q_med
sort race sex age
by race sex: egen totaldxi = total(dxi)
table race sex, c(min totaldxi)
	// Matches ELE main result

* dxi-weighted ex: first-pass, back-of-the-envelope average life lost 
* among those dying of police violence
* This measure understates life lost because it's based on ex including those deaths
* But more importantly, it overstates life lost because it uses ex for the 
* beginning of the age interval, when deaths occur throughout the interval
table race sex [iweight=dxi], c(mean ex_corrected)
	// A bit lower than Bui et al's (using different data)

* dxi-weighted median age: typical age of those dying in police encounters
table race sex [iweight=dxi], c(p50 age)
	// White and API: 35-39 category. Everyone else: 30-34.
	// The racial groups with the highest life expectancy being killed a bit older
		// may partially cancel out their otherwise greater life lost per death
	
* Comparison for Native and API populations: death-weighted alternative ex
* (not cause-deleted) based on e85 = Lx/dx
gen ex_alternative = 1/m if age==85
gen ex85_diff = ex_alternative - ex_corrected
table race sex if age==85, c(p50 ex_corrected p50 ex_alternative p50 ex85_diff)
	// The choice of e(85) estimation makes about 1-2 years per capita difference
gen ex85_popdiff = totaldxi * ex85_diff
table race sex if age==85, c(min ex85_popdiff)
	// The choice of e(85) estimation makes 111 years' difference in cohort life
		// lost (per 100,000) for Native men; 15 or under for others

* Multiple-decrement life table
ren q_med qxi
ren m_med mxi
gen qx = dx/lx
sum qx if age==85 // should be 1; is except for tiny rounding errors
replace qx = 1 if age==85

*** Associated single-decrement life table 
* (for single-decrement "not deaths to police force")
sort sex race age
gen R_negi = (dx - dxi)/dx
gen px_negi = (1-qx)^R_negi
replace px_negi = 0 if age==85

gen lx_negi = 100000 if age==0
replace lx_negi = lx_negi[_n-1] * px_negi[_n-1] if age==1
forvalues a=5(5)85 {
	replace lx_negi = lx_negi[_n-1] * px_negi[_n-1] if age==`a'
}

gen dx_negi = lx_negi * (1-px_negi)

gen n = 5 if age>=5
replace n = 1 if age==0
replace n = 4 if age==1

gen Lx_negi = n*lx_negi[_n+1] + dx_negi * nax if age<85
replace Lx_negi = lx_negi * ex_corrected / R_negi if age==85

gen Tx_negi = Lx_negi if age==85
forvalues a=80(-5)5 {
	replace Tx_negi = Lx_negi + Tx_negi[_n+1] if age==`a'
}
foreach a in 1 0 {
	replace Tx_negi = Lx_negi + Tx_negi[_n+1] if age==`a'
}
gen ex_negi = Tx_negi/lx_negi
list age lx_negi dx_negi Lx_negi Tx_negi ex_negi if race==1 & sex=="Male"

* true average life lost
gen lifelost_i_of_ind = ex_negi - nax
table race sex [iweight=dxi], c(mean lifelost_i_of_ind)

* Measure: total life-years lost per 100,000
* equals deaths per 100,000 * life lost (ex)
sort sex race age
gen lifelost_i_of = dxi * (ex_negi - nax)
gen cumlifelost_i_of = lifelost_i_of if age==0
replace cumlifelost_i_of = lifelost_i_of + cumlifelost_i_of[_n-1] if age==1
forvalues a=5(5)85 {
	replace cumlifelost_i_of = lifelost_i_of + cumlifelost_i_of[_n-1] if age==`a'
}
list age dxi ex_corrected lifelost_i_of cumlifelost_i_of if ///
	race==1 & sex=="Male" 
	// Cumulative life lost to officer force 3,772.1336 years for Black men
line lifelost_i_of age if race==1 & sex=="Male"
line cumlifelost_i_of age if race==1 & sex=="Male"

* save file with CDLT of officer force
save CDLT_of.dta, replace

* --------------------------------
* 2. Incorporating other fatal encounters
* --------------------------------

use fe_imp.dta, clear // includes 2013-2018 estimates 
	// and multiple imputations for force estimates (from ELE data)
* compare min and max
sort year sex race age
foreach var in officer_force vehicle suicide other {
	by year sex race age: egen tempmin`var' = min(`var')
	by year sex race age: egen tempmax`var' = max(`var')
}
sum temp*
drop temp*	
table age race sex, c(min officer_force max officer_force) 
	// range over years and imputations
	
* create total-race category
egen tempgroup = group(year sex age _imp)
expandby 2, by(tempgroup) gen(new)
replace race="Total" if new==1
drop tempgroup new
sort _imp year sex age race
foreach var in officer_force other suicide vehicle pop {
	replace `var' = . if race=="Total"
	by _imp year sex age: egen temp`var' = total(`var')
	replace `var' = temp`var' if race=="Total"
}
drop temp*

* make deaths and death rates -- averaging over years and imputations
* These average in that both deaths and populations are totaled over
	* all 10 imputations
gen alldeaths = officer_force + vehicle + suicide + other
gen forcevehicledeaths = officer_force + vehicle
gen mx_all = alldeaths/pop
gen mx_fv = forcevehicledeaths/pop
gen mx_of = officer_force / pop

* Visualize the mortality rates
gen bm = 1 if sex=="Male" & race=="African American"
scatter mx_all age if bm==1 || ///
	scatter mx_fv age if bm==1 || scatter mx_of age if bm==1

table race sex if age==25, c(p50 mx_all p50 mx_fv p50 mx_of)
table race sex if age==30, c(p50 mx_all p50 mx_fv p50 mx_of)

* Collapse over year and imputation
collapse (sum) deaths_all=alldeaths deaths_fv=forcevehicledeaths ///
	deaths_of=officer_force pop=pop, by(race sex age)
* Divide by the 10 imputations
* (Population is also divided because it also summed over imputations in the collapse)
foreach var in deaths_all deaths_fv deaths_of pop {
	replace `var' = `var'/10
}
	
gen bm = 1 if sex=="Male" & race=="African American"
foreach v in all fv of {
	gen mx_`v' = deaths_`v' / pop
	gen mx100k_`v' = mx_`v' * 100000
}
list age mx* if bm==1
table race sex, c(max mx100k_all max mx100k_fv max mx100k_of) // max death rates at any age

* Merge with main life table
	ren race racename
	gen race = 1 if racename=="African American"
	replace race = 2 if racename=="American Indian/AK Native"
	replace race = 3 if racename=="Asian/Pacific Islander"
	replace race = 4 if racename=="Latinx"
	replace race = 5 if racename=="White"
	replace race = 6 if racename=="Total"
	label define racelbl 1 "Black" 2 "Native" 3 "Asian/PI" ///
		4 "Latinx" 5 "White" 6 "Total"
	label values race racelbl
	
merge 1:1 sex race age using CDLT_of.dta, update
drop _merge	
sort sex race age
save CDLT_all_basefile.dta, replace

**** Cause-deleted life table for ALL police-associated deaths
drop qxi mxi *_negi totaldxi // drop most CDLT variables for excessive force only
					// These will be replaced by new variables for all fatal encounters
ren dxi dxi_of
** First, associated single-decrement table
* Make R_negi, based on dx_i, requires qx_i
gen mxi_all = mx_all
gen qxi_all = n * mxi_all / (1 + (n - nax)*mxi_all)
gen dxi_all = qxi_all * lx // true lx, not lx_i

* Table of total dx
sort race sex age
by race sex: egen totaldxi_all = total(dxi_all)
table race sex, c(min totaldxi_all) // black men: 141 per 100,000

* dxi-weighted ex: first-pass, back-of-the-envelope average life lost 
* among those dying of police violence
* This measure understates life lost because it's based on ex including those deaths
* But more importantly, it overstates life lost because it uses ex for the 
* beginning of the age interval, when deaths occur throughout the interval
table race sex [iweight=dxi_all], c(mean ex_corrected)

* dxi-weighted median age: typical age of those dying in police encounters
table race sex [iweight=dxi_all], c(p50 age)
	* As with officer use of force except Latinx men shift to 35-39 category

*** Associated single-decrement life table 
* (for single-decrement "not deaths to police excessive force")
sort sex race age
gen R_negi_all = (dx - dxi_all)/dx
gen px_negi_all = (1-qx)^R_negi_all
replace px_negi_all = 0 if age==85

gen lx_negi_all = 100000 if age==0
replace lx_negi_all = lx_negi_all[_n-1] * px_negi_all[_n-1] if age==1
forvalues a=5(5)85 {
	replace lx_negi_all = lx_negi_all[_n-1] * px_negi_all[_n-1] if age==`a'
}

gen dx_negi_all = lx_negi_all * (1-px_negi_all)

gen Lx_negi_all = n*lx_negi_all[_n+1] + dx_negi_all * nax if age<85
replace Lx_negi_all = lx_negi_all * ex_corrected / R_negi_all if age==85

gen Tx_negi_all = Lx_negi_all if age==85
forvalues a=80(-5)5 {
	replace Tx_negi_all = Lx_negi_all + Tx_negi_all[_n+1] if age==`a'
}
foreach a in 1 0 {
	replace Tx_negi_all = Lx_negi_all + Tx_negi_all[_n+1] if age==`a'
}
gen ex_negi_all = Tx_negi_all/lx_negi_all
list age lx_negi_all dx_negi_all Lx_negi_all Tx_negi_all ex_negi_all ///
	if race==1 & sex=="Male"

* true average life lost
gen lifelost_i_all_ind = ex_negi_all - nax
table race sex [iweight=dxi_all], c(mean lifelost_i_all_ind)	
	
* Measure: total life-years lost per 100,000
* equals deaths per 100,000 * life lost (ex)
sort sex race age
gen lifelost_i_all = dxi_all * lifelost_i_all_ind
	// estimate of life lost is based on life expectancy in the *absence* of these deaths
gen cumlifelost_i_all = lifelost_i_all if age==0
replace cumlifelost_i_all = lifelost_i_all + cumlifelost_i_all[_n-1] if age==1
forvalues a=5(5)85 {
	replace cumlifelost_i_all = lifelost_i_all + cumlifelost_i_all[_n-1] if age==`a'
}

* Examine life lost
list age dxi_all ex_corrected lifelost_i_all cumlifelost_i_all if bm==1
	// * For every 100,000 Black men, 5,696.3441 years of life lost
	
merge 1:1 sex race age using CDLT_of.dta, keepusing(lifelost* cumlife*)
drop _merge
sort sex race age	
line cumlifelost_i_of age if bm==1, legend(label(1 "Excessive force")) ///
		lwidth(thick) || ///
	line cumlifelost_i_all age if bm==1, ///
		legend(label(2 "All police-associated deaths")) ///
		lwidth(thick) ///
	xtitle("Age") ytitle("Cumulative life-years lost per 100,000") ///
	subtitle("Black Men's Life Years Lost to Police-Associated Deaths") ///
	saving(lifeyearslost_age.gph, replace)
	graph export lifeyearslost_age.png, replace
	
table race sex, c(max cumlifelost_i_of max cumlifelost_i_all)

* Proportion of life lost to police encounters that is due to officer force
gen prop_of = cumlifelost_i_of / cumlifelost_i_all if age==85
table race sex, c(min prop_of)

save CDLT_all.dta, replace

* --------------------------------
* 3. Bar graphs
* --------------------------------
* Check population sizes and make race variable ordered by size
bysort race sex: egen totalpop = total(pop)
table race sex, c(min totalpop)
	* Smallest to largest: Native, Asian/PI, Black, Latinx, white, total
ren race raceorig
recode raceorig (2=1) (3=2) (1=3), gen(race)
label define racelbl_final 1 "Native" 2 "Asian/PI" 3 "Black" ///
		4 "Latinx" 5 "White" 6 "Total"
label values race racelbl_final
tab race raceorig
	
* Set up variables
gen cumall = cumlifelost_i_all if age==85
gen cumof = cumlifelost_i_of if age==85
clonevar cumallrev = cumall
replace cumallrev = -cumall if sex=="Female"
clonevar cumofrev = cumof
replace cumofrev = -cumof if sex=="Female"

* Population pyramid-style graph
twoway bar cumallrev race if sex=="Male", horizontal ///
		color("204 0 51") ///
		/// blabel(race, position(inside)) ///
		legend(label(1 "All police-associated, Men")) || ///
	bar cumallrev race if sex=="Female", horizontal ///
		color("0 51 153") ///
		legend(label(2 "All police-associated, Women")) || ///
	bar cumofrev race if sex=="Male", horizontal ///
		color("204 102 51") ///
		/// blabel(race, position(inside)) ///
		legend(label(3 "Officer force only, Men")) || ///
	bar cumofrev race if sex=="Female", horizontal ///
		color("0 153 153") ///
		legend(label(4 "Officer force only, Women")) ///
	xlabel(-600 "600" ///
		0(1000)6000) xtick(-600 -300 0(1000)6000) ///
	ytitle("") ///
	title("Life lost, per 100,000 people, to police encounters") ///
	legend(on) ///
	saving(lifelost_bar.gph, replace)
	graph export lifelost_bar.png, replace

	/*	
graph hbar cumall race if sex=="Male",  ///
		/// color("204 0 51") ///
		blabel(race, position(inside)) ///
		legend(label(1 "All police-associated, Men")) || ///
	graph hbar cumall race if sex=="Female",  ///
		/// color("0 51 153") ///
		legend(label(2 "All police-associated, Women")) || ///
	graph hbar cumof race if sex=="Male",  ///
		/// color("204 102 51") ///
		/// blabel(race, position(inside)) ///
		legend(label(3 "Excessive force only, Men")) || ///
	graph hbar cumof race if sex=="Female",  ///
		/// color("0 153 153") ///
		legend(label(4 "Excessive force only, Women")) ///
	xlabel(-300 "300" 0(1000)6000) xtick(-300 0(1000)6000) ///
	ytitle("") ///
	bar(1, color("204,0,51")) ///
	bar(3, color("204,102,51")) ///
	bar(2, color("0,51,153")) ///
	bar(4, color("0,153,153")) ///
	title("Life lost, per 100,000 people, to police encounters") ///
	saving(lifelost_bar2.gph, replace)
	graph export lifelost_bar2.png, replace	
	
	
graph hbar cumof cumall, over(race, reverse) over(sex) ///
	legend(label(1 "Excessive force only")) ///
	legend(label(2 "All police-associated")) ///
	title("Life lost, per 100,000 people, to police encounters")

graph hbar cumof cumall if sex=="Male", over(race, reverse) ///
	legend(label(1 "Excessive force only")) ///
	legend(label(2 "All police-associated")) ///
	title("Life lost, per 100,000 men, to police encounters")
	*/
	
* Separate panels with full labels and legends:	
graph hbar cumall cumof if sex=="Male", over(race, /* reverse */ ) ///
	blabel(position(inside), color(white)) ///
	legend(label(2 "Officer force only")) ///
	legend(label(1 "All police-associated deaths")) ///
	legend(on order(2 1) size(medsmall)) ///
	title("Men's life-years lost, per 100,000 men, to police encounters", ///
		size(medlarge))	///
	saving(lifelost_hbar_men.gph, replace)
	graph export lifelost_hbar_men.png, replace
	
graph hbar cumall cumof if sex=="Male", ///
	over(race, lab(nolab) /* reverse */ ) ///
	blabel(position(inside), color(white)) ///
	legend(label(2 "Officer force only")) ///
	legend(label(1 "All police-associated deaths")) ///
	legend(on order(2 1) size(medsmall)) ///
	title("Men", ///
		size(medlarge))	///
	saving(lifelost_hbar_men_noxshorttitle.gph, replace)
	
graph hbar cumall cumof if sex=="Male", ///
	over(race, /* reverse */ ) ///
	blabel(position(inside), color(white)) ///
	legend(label(2 "Officer force only")) ///
	legend(label(1 "All police-associated deaths")) ///
	legend(on order(2 1) size(medsmall)) ///
	title("Men", ///
		size(medlarge))	///
	saving(lifelost_hbar_men_shorttitle.gph, replace)	
	
graph hbar cumall cumof if sex=="Female", over(race, /* reverse */) ///
	blabel(position(inside), color(white)) ///
	legend(label(2 "Officer force only")) ///
	legend(label(1 "All police-associated deaths")) ///
	legend(/* on */ order(2 1) size(medsmall)) ///
	title("Women's life-years lost, per 100,000 women, to police encounters", ///
		size(medlarge))	///
	saving(lifelost_hbar_women.gph, replace)
	graph export lifelost_hbar_women.png, replace	
	
graph hbar cumall cumof if sex=="Female", over(race, /* reverse */) ///
	blabel(position(inside), color(white)) ///
	legend(label(2 "Officer force only")) ///
	legend(label(1 "All police-associated deaths")) ///
	legend(/* on */ order(2 1) size(medsmall)) ///
	title("Women", ///
		size(medlarge))	///
	saving(lifelost_hbar_women_shorttitle.gph, replace)

grc1leg lifelost_hbar_men_shorttitle.gph lifelost_hbar_women_shorttitle.gph, ///
	title("Life Years Lost to Police Encounters per 100,000 People") ///
	note("Note that the scale for men is 10 times greater than for women", ///
		margin(l+10.5)) ///
	saving(lifelost_hbar_menwomen.gph, replace)
	graph export lifelost_hbar_menwomen.png, replace
	graph export lifelost_hbar_menwomen.tif, width(3500) replace
		
* Panels with women on reversed axis and labels, assembled into a 
* population pyramid-style graph:
		
graph hbar cumall cumof if sex=="Male", ///
	over(race, lab(nolab) /* reverse */ ) ///
	blabel(position(inside), color(white)) ///
	legend(label(2 "Officer force only")) ///
	legend(label(1 "All police-associated deaths")) ///
	legend(off order(2 1) size(medsmall)) ///
	title("Men's life-years lost, per 100,000 men, to police encounters", ///
		size(medlarge))	///
	saving(lifelost_hbar_men_nox.gph, replace)

graph hbar cumallrev cumofrev if sex=="Female", over(race, /* reverse */) ///
	blabel(position(inside), color(white)) ///
	xalternate ///
	legend(label(2 "Officer force only")) ///
	legend(label(1 "All police-associated deaths")) ///
	legend(/* on */ order(2 1) size(medsmall)) ///
	ylabel(0 "0" -200 "200" -400 "400" -600 "600") ///
	title("Women", size(medlarge))	///
	saving(lifelost_hbar_women_yalt.gph, replace)
	
graph hbar cumallrev cumofrev if sex=="Female", ///
	over(race, lab(nolab) /* reverse */) ///
	blabel(position(inside), color(white)) ///
	xalternate ///
	legend(label(2 "Officer force only")) ///
	legend(label(1 "All police-associated deaths")) ///
	legend(off order(2 1) size(medsmall)) ///
	ylabel(0 "0" -200 "200" -400 "400" -600 "600") ///
	title("Women", ///
		size(medlarge))	///
	saving(lifelost_hbar_women_yalt_nox.gph, replace)	
		
grc1leg lifelost_hbar_women_yalt.gph lifelost_hbar_men_noxshorttitle.gph, ///
	title("Life Years Lost per 100,000 People") ///
	note("Note that the scale for men is 10 times greater than for women") ///
	saving(lifelost_hbarpyramid_menwomen.gph, replace)
	graph export lifelost_hbarpyramid_menwomen.png, replace
		
* --------------------------------
* 4. Lifetime deaths, applying these rates to the 2017 birth cohort
* --------------------------------
* Black male births:
* black births in 2017 estimated at 560,715 by Kaiser
* from https://www.kff.org/other/state-indicator/births-by-raceethnicity/
* black sex ratio in 2018 estimated at 102.9 (0.507)
* from https://www.cdc.gov/nchs/data/nvsr/nvsr68/nvsr68_13-508.pdf
scalar bm_births = 560715 * (1029/2029) // Black births * black SRB in proportion form
scalar bm_popmult = bm_births / 100000
foreach var in of all {
	sum cum`var' if bm==1
	scalar bm_cum`var' = r(mean)
	scalar bm_poplifelost_`var' = bm_cum`var' * bm_popmult
}
scalar list

* Make cumulative deaths over age
sort sex race age
foreach var in all of {
	gen cumdxi_`var' = dxi_`var' if age==0
	replace cumdxi_`var' = dxi_`var' + cumdxi_`var'[_n-1] if age==1
	forvalues a=5(5)85 {
		replace cumdxi_`var' = dxi_`var' + cumdxi_`var'[_n-1] if age==`a'
	}
	* Use those to estimate total deaths by age for the 2017 birth cohort
	gen cumdeaths_pop_`var' = cumdxi_`var' * bm_popmult if bm==1
}
list age dxi* cumdxi* cumdeaths* if bm==1

list cumall cumof cumdeaths* if bm==1 & age==85

save CDLT_all.dta, replace

log close
* done!
