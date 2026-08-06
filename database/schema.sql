CREATE TABLE compound (
    inchikey      VARCHAR(27) PRIMARY KEY,
    smiles        TEXT,
    name          VARCHAR(255)
);


CREATE TABLE chembl (
    chembl_id     VARCHAR(20) PRIMARY KEY,
    inchikey      VARCHAR(27) NOT NULL,
    CONSTRAINT fk_chembl_compound
        FOREIGN KEY (inchikey) REFERENCES compound (inchikey)
);


CREATE TABLE uniprot (
    uniprot_id    VARCHAR(20) PRIMARY KEY,
    entrez_gene   VARCHAR(50),
    hgnc          VARCHAR(50),
    species       VARCHAR(100)
);


CREATE TABLE target (
    target_id     SERIAL PRIMARY KEY,
    type          VARCHAR(50) NOT NULL,
    name          VARCHAR(255),
    CONSTRAINT ck_target_type
        CHECK (type IN ('protein', 'complex', 'ppi', 'family'))
);


CREATE TABLE target_uniprot (
    target_id     INTEGER NOT NULL,
    uniprot_id    VARCHAR(20) NOT NULL,
    CONSTRAINT pk_target_uniprot PRIMARY KEY (target_id, uniprot_id),
    CONSTRAINT fk_tu_target
        FOREIGN KEY (target_id) REFERENCES target (target_id),
    CONSTRAINT fk_tu_uniprot
        FOREIGN KEY (uniprot_id) REFERENCES uniprot (uniprot_id)
);

-- source_db is the resource, source is the record inside it (a paper, an
-- internal report, a release). xref_id is the prefix that turns the
-- bioactivity.source_xref of a single measurement into a link, so for a
-- literature source it is https://doi.org/ and source_xref is the DOI.
-- leave xref_id empty if the source has nothing resolvable.
CREATE TABLE bioactivity_source (
    source_id     SERIAL PRIMARY KEY,
    source_db     VARCHAR(255) NOT NULL,
    source        VARCHAR(255),
    xref_id       VARCHAR(255),
    CONSTRAINT uq_bioactivity_source UNIQUE (source_db, source)
);


CREATE TABLE bioactivity_group (
    id            SERIAL PRIMARY KEY,
    inchikey      VARCHAR(27) NOT NULL,
    target_id     INTEGER NOT NULL,
    moa           VARCHAR(255) NOT NULL DEFAULT '',
    CONSTRAINT fk_bag_compound
        FOREIGN KEY (inchikey) REFERENCES compound (inchikey),
    CONSTRAINT fk_bag_target
        FOREIGN KEY (target_id) REFERENCES target (target_id),
    CONSTRAINT uq_bioactivity_group
        UNIQUE (inchikey, target_id, moa)
);


CREATE TABLE bioactivity (
    id                SERIAL PRIMARY KEY,
    inchikey          VARCHAR(27) NOT NULL,
    target_id         INTEGER NOT NULL,
    moa               VARCHAR(255) NOT NULL DEFAULT '',
    bioactivity_type  VARCHAR(50), 
    relation          VARCHAR(5),
    value             NUMERIC,
    unit              VARCHAR(50),
    assay_type        VARCHAR(50),
    assay_description VARCHAR(255),
    cell_line         VARCHAR(100),
    concentration     NUMERIC,
    concentration_unit VARCHAR(50),
    source_id         INTEGER,
    source_xref       VARCHAR(100),
    CONSTRAINT ck_relation
        CHECK (relation IS NULL OR relation IN ('=', '>', '<', '>=', '<=', '~')),
    CONSTRAINT fk_ba_compound
        FOREIGN KEY (inchikey) REFERENCES compound (inchikey),
    CONSTRAINT fk_ba_target
        FOREIGN KEY (target_id) REFERENCES target (target_id),
    CONSTRAINT fk_ba_source
        FOREIGN KEY (source_id) REFERENCES bioactivity_source (source_id),
    CONSTRAINT fk_ba_group
        FOREIGN KEY (inchikey, target_id, moa)
        REFERENCES bioactivity_group (inchikey, target_id, moa)
);

CREATE INDEX idx_target_uniprot_acc    ON target_uniprot (uniprot_id);
CREATE INDEX idx_bioactivity_compound  ON bioactivity (inchikey);
CREATE INDEX idx_bioactivity_target    ON bioactivity (target_id);
CREATE INDEX idx_bioactivity_source    ON bioactivity (source_id);
CREATE INDEX idx_bioactivity_scope     ON bioactivity (bioactivity_type, unit);
CREATE INDEX idx_compound_name         ON compound (name);
CREATE INDEX idx_chembl_inchikey       ON chembl (inchikey);


-- target and uniprot have no column in common, the link is target_uniprot.
-- this view does both joins so you get one flat frame instead of three.
-- LEFT so a target with no accession yet still shows up.
CREATE VIEW target_flat AS
SELECT t.target_id, t.type, t.name,
       tu.uniprot_id, u.hgnc, u.species, u.entrez_gene
  FROM target t
  LEFT JOIN target_uniprot tu ON tu.target_id = t.target_id
  LEFT JOIN uniprot u ON u.uniprot_id = tu.uniprot_id;



-- ---------------------------------------------------------------------------
-- Added for sources that carry more than a measurement (chemicalprobes.org/).
--
-- Everything above this line is unchanged. Nothing here is specific to one
-- source: every table names the resource that made the claim, so two sources
-- may disagree and both survive, which is the rule bioactivity already uses.
-- ---------------------------------------------------------------------------

-- what a source concludes about a compound as a probe. one row per compound
-- per source, so opnMe recommending what the portal rejected is two rows and
-- not a conflict. the flags are computable from the structure rather than
-- opinions, but a source states them, so they are stored as stated.
CREATE TABLE probe_assessment (
    inchikey           VARCHAR(27) NOT NULL,
    source_id          INTEGER NOT NULL,
    verdict            VARCHAR(20),
    pains              VARCHAR(3),
    toxicophore        VARCHAR(3),
    rating_in_cell     NUMERIC,
    rating_in_organism NUMERIC,
    rating_count       INTEGER,
    published_date     VARCHAR(20),
    CONSTRAINT pk_probe_assessment PRIMARY KEY (inchikey, source_id),
    CONSTRAINT ck_verdict     CHECK (verdict IS NULL OR verdict IN ('recommended', 'unsuitable')),
    CONSTRAINT ck_pains       CHECK (pains IS NULL OR pains IN ('Yes', 'No')),
    CONSTRAINT ck_toxicophore CHECK (toxicophore IS NULL OR toxicophore IN ('Yes', 'No')),
    CONSTRAINT fk_pa_compound FOREIGN KEY (inchikey) REFERENCES compound (inchikey),
    CONSTRAINT fk_pa_source   FOREIGN KEY (source_id) REFERENCES bioactivity_source (source_id)
);


-- an identifier for a compound in some other resource. chembl has its own table
-- because it came first; this is the same idea for the rest, so the next
-- resource is one value in a CHECK and not a migration.
CREATE TABLE compound_xref (
    resource      VARCHAR(50) NOT NULL,
    xref          VARCHAR(100) NOT NULL,
    inchikey      VARCHAR(27) NOT NULL,
    CONSTRAINT pk_compound_xref PRIMARY KEY (resource, xref),
    CONSTRAINT ck_resource CHECK (resource IN ('canSAR', 'PubChem', 'DrugBank', 'ZINC')),
    CONSTRAINT fk_cx_compound FOREIGN KEY (inchikey) REFERENCES compound (inchikey)
);


-- where a source puts a target in its own taxonomy. two levels, and both are
-- multi-valued: one target is filed under several classes by different records,
-- so ordinal keeps them all rather than letting the last one win.
CREATE TABLE target_class (
    target_id     INTEGER NOT NULL,
    source_id     INTEGER NOT NULL,
    level         VARCHAR(10) NOT NULL,
    ordinal       INTEGER NOT NULL DEFAULT 0,
    value         VARCHAR(255) NOT NULL,
    CONSTRAINT pk_target_class PRIMARY KEY (target_id, source_id, level, ordinal),
    CONSTRAINT ck_level CHECK (level IN ('class', 'subclass')),
    CONSTRAINT fk_tc_target FOREIGN KEY (target_id) REFERENCES target (target_id),
    CONSTRAINT fk_tc_source FOREIGN KEY (source_id) REFERENCES bioactivity_source (source_id)
);


-- a dose given to an animal. not a potency: no endpoint, nothing to compare,
-- and no target, so it cannot be a bioactivity row. one row per dose, and
-- dose_raw keeps the string it was read out of.
CREATE TABLE in_vivo_dose (
    id            SERIAL PRIMARY KEY,
    inchikey      VARCHAR(27) NOT NULL,
    source_id     INTEGER NOT NULL,
    organism      VARCHAR(100),
    dose_value    NUMERIC,
    dose_unit     VARCHAR(50),
    route         VARCHAR(20),
    dose_raw      VARCHAR(255),
    CONSTRAINT ck_route CHECK (route IS NULL OR route IN ('IV', 'PO', 'IP', 'SC', 'topical')),
    CONSTRAINT fk_ivd_compound FOREIGN KEY (inchikey) REFERENCES compound (inchikey),
    CONSTRAINT fk_ivd_source   FOREIGN KEY (source_id) REFERENCES bioactivity_source (source_id)
);


-- the literature a source cites for a compound, where it cites a reading list
-- rather than a paper per measurement. raw is kept because 88% of these urls
-- cannot be rebuilt from xref_id and source_xref.
CREATE TABLE compound_reference (
    id            SERIAL PRIMARY KEY,
    inchikey      VARCHAR(27) NOT NULL,
    source_id     INTEGER NOT NULL,
    xref_id       VARCHAR(255),
    source_xref   VARCHAR(400),
    raw           TEXT,
    CONSTRAINT fk_cr_compound FOREIGN KEY (inchikey) REFERENCES compound (inchikey),
    CONSTRAINT fk_cr_source   FOREIGN KEY (source_id) REFERENCES bioactivity_source (source_id)
);


-- the escape hatch, for a fact a source states that has no column anywhere.
-- one property per row, ordinal so a property can hold a list. keep it small:
-- anything every source has deserves a column instead.
CREATE TABLE compound_annotation (
    inchikey      VARCHAR(27) NOT NULL,
    source_id     INTEGER NOT NULL,
    property      VARCHAR(100) NOT NULL,
    ordinal       INTEGER NOT NULL DEFAULT 0,
    value         TEXT,
    CONSTRAINT pk_compound_annotation PRIMARY KEY (inchikey, source_id, property, ordinal),
    CONSTRAINT fk_ca_compound FOREIGN KEY (inchikey) REFERENCES compound (inchikey),
    CONSTRAINT fk_ca_source   FOREIGN KEY (source_id) REFERENCES bioactivity_source (source_id)
);


-- one curator worklist: a record a source published that did not become a row.
-- a number that could be read but not trusted, or a compound with no structure
-- at all. stage says which, reason says why, raw keeps what it came from.
CREATE TABLE rejected_record (
    id            SERIAL PRIMARY KEY,
    source_id     INTEGER NOT NULL,
    stage         VARCHAR(20) NOT NULL,
    reason        VARCHAR(255) NOT NULL,
    label         VARCHAR(255),
    inchikey      VARCHAR(27),
    target_id     INTEGER,
    raw           TEXT,
    CONSTRAINT ck_stage CHECK (stage IN ('compound', 'bioactivity')),
    CONSTRAINT fk_rr_source   FOREIGN KEY (source_id) REFERENCES bioactivity_source (source_id),
    CONSTRAINT fk_rr_compound FOREIGN KEY (inchikey) REFERENCES compound (inchikey),
    CONSTRAINT fk_rr_target   FOREIGN KEY (target_id) REFERENCES target (target_id)
);

CREATE INDEX idx_probe_assessment_verdict ON probe_assessment (verdict);
CREATE INDEX idx_compound_xref_inchikey   ON compound_xref (inchikey);
CREATE INDEX idx_target_class_target      ON target_class (target_id);
CREATE INDEX idx_in_vivo_dose_compound    ON in_vivo_dose (inchikey);
CREATE INDEX idx_compound_reference_cmpd  ON compound_reference (inchikey);
CREATE INDEX idx_rejected_record_source   ON rejected_record (source_id);


-- one row per compound per source: the assessment with its identifiers, which
-- is the sheet a curator actually wants. no chembl join, because two compounds
-- carry two ChEMBL ids and it would fan out.
CREATE VIEW probe_flat AS
SELECT a.inchikey, c.name, c.smiles, s.source_db, a.verdict, a.pains, a.toxicophore,
       a.rating_in_cell, a.rating_in_organism, a.rating_count, a.published_date,
       (SELECT x.xref FROM compound_xref x
         WHERE x.inchikey = a.inchikey AND x.resource = 'canSAR') AS cansar_id
  FROM probe_assessment a
  JOIN compound c ON c.inchikey = a.inchikey
  JOIN bioactivity_source s ON s.source_id = a.source_id;


-- target_class is two rows deep for a reason, but "what class is BRD4" should
-- not need three joins to answer.
CREATE VIEW target_class_flat AS
SELECT t.target_id, t.name, u.hgnc, tc.level, tc.value, s.source_db
  FROM target_class tc
  JOIN target t ON t.target_id = tc.target_id
  JOIN bioactivity_source s ON s.source_id = tc.source_id
  LEFT JOIN target_uniprot tu ON tu.target_id = t.target_id
  LEFT JOIN uniprot u ON u.uniprot_id = tu.uniprot_id;
