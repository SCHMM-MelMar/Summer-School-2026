# General Discussion

## 1. User needs
Everybody has read the user needs. Following point five (specifying user needs that users come up with themselves) two people present use cases that could be unknown to others:

- The first requests the database to contain type of affinity
- The second requires thought about notation of phenotype in the database

## 2. Database Schemes
### Db1:
Compound, target and source tables exist in parallel and have IDc, IDp and IDs ids respectively.
A combination table of these exists with IDc, IDp and IDs specified, with the interaction information noted in this combination table.

Because targets can be associated to different classes (PROTAC, complex and conventional protein binding differ fundamentally), each class gets an information table that is related one-to-one to entries in the target table. These information tables are to reduce redundancy in the target table.

### Db2:
Central in the database is a bioactivity table that combines compound and target and their associated properties (MoA, affinity etc.). This also contains the source of the bioactivities.

This database is linked to a compound table that shows all properties that are compound-specific. Furthermore, this is linked to a target table that shows all properties that are associated to the target.

The target table can contain multiple uniprot entries, these are shown in a separate uniprot table, where a one-to-many relationship associates the target and uniprot tables.

All bio-activities are further associated with a bioactivity-group table, which is necessary to ensure that multiple entries in the bio-activities table that have the same compound-id, target-id and method of action, are able to exist as separate entries in the database.

The compound table is related to a table that associates compound-ids to ids in other databases. Normally, this could be avoided, and one would only use chembl-id. Unfortunately, chembl-id’s only exist for publicly available compounds. Therefore creating a table that shows the source and is able to associate to non-public databases (such as Roche’s own identifiers) is of importance.

Noteworthy is that expert knowledge is missing in this design. Suggestion is that expert intuition would be added in another table which allows usernames and comments to be associated to both targets as well as to compounds.

## To Be Continued
GROUP DISCUSSION POSTPONED TILL AFTER LUNCH AND TOUR
After lunch reconvene to discuss the following:
- Scheme discussion and determination
- Plan determination
- Group task assignment
