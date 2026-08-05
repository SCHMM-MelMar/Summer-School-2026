# General Discussion

## 1. User needs
Everybody has read the user needs. Following point five (specifying user needs that users come up with themselves) two people present use cases that could be unknown to others:

- The first requests the database to contain type of affinity
- The second requires thought about notation of phenotype in the database

## 2. Database Schemes
### Db1:

![Fig 1: database scheme for db1](/discussions/images/db1.png)

Compound, target and source tables exist in parallel and have IDc, IDp and IDs ids respectively.
A combination table of these exists with IDc, IDp and IDs specified, with the interaction information noted in this combination table.

Because targets can be associated to different classes (PROTAC, complex and conventional protein binding differ fundamentally), each class gets an information table that is related one-to-one to entries in the target table. These information tables are to reduce redundancy in the target table.

### Db2:

![Fig 2: database scheme for db2](/discussions/images/db2.png)

Central in the database is a bioactivity table that combines compound and target and their associated properties (MoA, affinity etc.). This also contains the source of the bioactivities.

This database is linked to a compound table that shows all properties that are compound-specific. Furthermore, this is linked to a target table that shows all properties that are associated to the target.

The target table can contain multiple uniprot entries, these are shown in a separate uniprot table, where a one-to-many relationship associates the target and uniprot tables.

All bio-activities are further associated with a bioactivity-group table, which is necessary to ensure that multiple entries in the bio-activities table that have the same compound-id, target-id and method of action, are able to exist as separate entries in the database.

The compound table is related to a table that associates compound-ids to ids in other databases. Normally, this could be avoided, and one would only use chembl-id. Unfortunately, chembl-id’s only exist for publicly available compounds. Therefore creating a table that shows the source and is able to associate to non-public databases (such as Roche’s own identifiers) is of importance.

Noteworthy is that expert knowledge is missing in this design. Suggestion is that expert intuition would be added in another table which allows usernames and comments to be associated to both targets as well as to compounds.

## Database Architecture Decision
P&D database is somewhat similar to the tables as proposed in database 2. The full database is 70 tables in size. It includes information on pathways and ontologies, and cross-references. “The basic structure is similar, the functionality and the data is there. It is more complex than we want it” -Dennis. Dennis recommends that we either create our own database and populate it with the databases that we have looked at ourselves, or we strip down the database of P&D.

David does not want to discus inclusion of phenotype for simplicity purposes, as the explanation of what phenotype means is complicated and non-structured data for in the database. Furthermore, it is unlikely that phenotype is important for the user cases and analysis.

Two options have been presented. Either the P&D dataset would be simplified, as their scaffold is largely similar to the targeted database, or a self-designed database is to be built based on the described architecture in db2 (and in part in db1).

Danya shows being a proponent for building both one from scratch and comparing it to P&D, to get an idea if our database has a usage benefit. This would be a QC step at the end. Nicky states that she agrees with the sentiment, especially as building it ourselves pushes for more understanding and learning. Lukas thinks this is in line with the intent of the Summer School, which might differ from the intended metrics. His reasoning found no further discourse (not necessarily agreement, but no outspoken opponents). 

11 people stated to have a preference to the self-built database, thus was decided to continue by building it ourselves.

Continuation of the requirements of the database and its columns was not recorded in minutes.

## Task Separation
Scheme Design (requires a representative of each database):

Database Masters:
ChemicalProbes.org: Nicky, Ben
SPARK: Evegenie, Daniel
MoA: Dennis, Andrea
P&D: Raquel, Brennan
OpmMe: Danya
SGC: Lukas
ReFrame: Andreea, Maedeh

IT Team:
Dennis, Danya, Melanie, Nina, Leander


## Day Goal
Prototype Database
Schedule Preprocessing Timelines

Will reconvene at the end of the day to discuss the architecture.
