== Change Log ==


----

## config.15.covid.mda

Main changes:
-	Mosquito seasonality annual curve
-	Reduction in treatment rate
-	MDA age limit


TODO:

-   y. treatment reduction 100 to 50 % in May

-   MDA on < 18 age

-   treatment.malaria_post, check 3 rates, use latest

-   y. check drug values are same as MDA paper

-   y. read old agestruct

-   y. use 1500 village regions file



### data

-	deathprob_v5_bmgf.txt
-	regions_newdata_mmc.txt, 1500 villages
-	agestruct15.txt

-	mosq season types


### config

#### structure modifications

-   Based on config.13
-   <s> mosquito.seasonality.data_file.annual_repeat </s>
-   mosquito.seasonality.data_file.horizontal_shift_days, represent the delay between a curve showing mosquito season and a curve showing new cases. About 40 days (1. incubation ~ 14 days, 2. liver to blood ~ 5 days, 3. blood to infectious (case) ~ 15 days, 4. a few days on feeding cycle)

-	Added treatment.mda_max_age
-	Added treatment.mda_min_age

-   removed human.itn_probability

-	Added fixed_routes option to MDAs


#### values

-   mosquito.seasonality.use_data_file = true
-   mosquito.seasonality.data_file.name = four types
-   mosquito.seasonality.data_file.horizontal_shift_days = 40


-   intervention.itn.init_probability_per_human = 0.75
-   intervention.itn. = 0.95 * 0.95 * 0.5

-   treatment.malaria_post.establish_malaria_post_on_day = 0
-   treatment.treatment_rate_mda = 0.5333

-   infection.func_m2h.version = 2

-   phi_a = 0.2

-   func_clinical_probability_blood_to_infectious.version = 4

-   mellow = 10 days

-   drug.allow_act_to_a = false
-   drug.possible_to_have_drug_a_only = false

-	attractiveness on
-	func_update_susceptibility


### Code

**New**

-   util::read_vector_from_file
-   VM::func_beta_to_beta_season_given_curve

-	


**Updated**

-   VM::step_update_beta_season_of

-	VM::treatment_by_malaria_post,  add reduction

-	VM::treatment_by_mda, add age limits



**Fixed**

-	BS::init_from_file