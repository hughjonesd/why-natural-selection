
Notes for census 2011 weightings
================================

* Data comes from https://www.ons.gov.uk/census/2011census/2011censusdata


Ideal breakdown
===============

Within ethnicity of "White British":
- i.e. NATIDPBUK (but do we include people with dual self-ascribed national 
  identity)

sex, age, education, household type, marital status, tenure at individual level. 
- n_in_household could match to SIZEHUK11
- f.6141 counts up to 4 others in the same household, and their relationships
  from that we could derive MARSTAT or LARPUK11 equivalents, as well as
  DPCFAMUK11 (dependent children in family)
- f.670 (type of accommodation) and f.680 (own or rent) are analogues to
  TENHUK11 and TYPEACCOM
- education will be complex... :-)
- marital status: MARSTAT (doesn't include cohabitation); LARPUK11 (cohabitation
  has more priority)

We could also match on geographies. That still leaves issues, because so many
areas are "left out". But it would weight us correctly within the geographies
that we have.

Since many of our couples have non-dependent children, it probably makes sense
to use census variables comparing simply living with partner/not, rather than
using presence or absence of children in the cells.

Household type makes a *huge* difference in the GHS weighting, so we probably
want to confirm this via LARPUK in the census data.

# Existing tables that would be useful:

DC1101EW - marital status by sex by age, England and Wales
DC1102EW - living arrangements ("reference person") by sex by age, E&W 
  - this is the same as DC1108EW but the sample is "all household reference persons"
  - we'd probably want to use "all persons" as we aren't sampling by household
DC1104EW - residence type by sex by age
DC1108EW - living arrangements by sex by age, E&W (but down to MSOA level, with
  some privacy protection swapping built in)
DC1109EW - household composition by sex by age, E&W (household composition
  is roughly married/cohabiting/lone parent/single x presence of dependent 
  children)
DC1201EW - household composition by ethnic group of "household reference person"
DC2101EW - ethnic group by sex by age *probably the simplest for MSOA level*
  - nb there is also LC2101EW where the sample is "all usual residents"
DC2102EW - national identity by sex by age
DC2103EW - country of birth by sex by age
DC4101EW - tenure by household composition
DC4201EW - tenure by ethnic group by age
DC5102EW - highest level of qualification by sex by age
DC5202EW - hlq by ethnic group by age
DC5801EW - as above, but also with year of arrival in the UK (so could do UK-born)

# What we can't do: 

- age at first live birth except when there are dependent
children in the family. (So, only for children under 16, since there's selection
on the 16-18 year olds. So, only among people having children after about 1990.)

- Number of siblings, parents etc.

- Income, directly.

- We can't differentiate cohabiting from married, so may as well use 
  "living arrangements" not "marital status"

Variables
=========

sex
* SEX (1 male, 2 female)

age
* AGE

education
* QUALS is the standard variable
* Derived results are:
  - QUALS01 (1-4 GCSEs/GCE/O level)
  - QUALS02 (NVQ level 1 or GNVQ)
  - QUALS03 (5+ GCSEs/GCE/O level)
  etc.  

ethnicity (closest equivalent to genetic_ethnic_grouping??? or to self-identification???)
* COB - country of birth - England/Scotland/Wales/NI/Rep of I/"elsewhere" + writein
  - Probably use
* PSPTEL - passports held - one or more of UK/Irish/"other" + writein/None
* NATIDPBUK - a self-assessment "how would you describe your national identity"
  - British, British and other, No British identity


age at first live birth

number of children
DPCHUK11 - dependent children in household (combines age and number)
DPCFAMUK11 - dependent children in *family* i.e. of a couple, lone parent, or 
  grandparent(s)

number of siblings?

income
NSSHUK11 - socioeconomic position. Roughly, class A-E, but with more detail.
OCCPUK113 - ditto, more broad, relates to SOC2010 highest level of hierarchy.
OCCPUK112 - more detail on same
OCCPUK111 - yet more
* Note: UKBB data is based on "SOC2000". At the highest level, the categories
  are very similar between SOC2010 and SOC2000: see 
  https://www.ons.gov.uk/methodology/classificationsandstandards/standardoccupationalclassificationsoc/soc2010/soc2010volume1structureanddescriptionsofunitgroups
  - in fact the 9 "major groups" are the same and "submajor groups" look very
    similar

household tenure
TENHUK11 - owned outright, with a mortgage, rented private/social...

household type
MARSTAT - marital status on census day. Single/married/separated/divorced/widowed
plus analogues from a civil partnership.
SIZHUK11 - number of people in household.
AHTHUK11 - married, cohabiting, lone parent etc.
LARPUK11 - living arrangements. Here  "cohabiting takes priority over other
  categories. For example, if a person is divorced and cohabiting, then in results
  for living arrangements they are classified as cohabiting."
ADTHUK11 - number of adults in household. (18+, or 16-18 but not in education - i.e.
  not a "dependent child")
P17HUK11 - number of people over 17 in the household.




Other possibles
===============

BEDROOMS   - number of bedrooms (bedsits/studios = 1)
TYPACCOM   - detached/semi/terrace/flat etc.
CARSNO     - cars or vans in household (0-10)
HOURS      - hours worked
MAINLANG   - first/preferred language
ACTLW - last week were they working, unemployed, retired, student etc.
RESIDENCE_TYPE = C (communal) or H (household)
INDUSTRY - where someone works, 22 categories from SIC2007
INDGEPUK11 - goes deeper into the same
MIGOAPW11 - has a person moved within the past year, and how far?
