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


-- the portal publishes compounds it has judged unsuitable as probes, with no
-- target and no measurement, so the only thing that explains them is this
-- verdict. name, smiles and the ChEMBL ids are not repeated here: inchikey is a
-- foreign key onto compound, which already holds them. source_id says who is
-- making the judgement, the same way a measurement does.
CREATE TABLE unsuitable (
    inchikey           VARCHAR(27) PRIMARY KEY,
    portal_path        VARCHAR(255),
    published_date     VARCHAR(20),
    pains              VARCHAR(3),
    toxicophore        VARCHAR(3),
    cansar_id          VARCHAR(20),
    rating_in_cell     NUMERIC,
    rating_in_organism NUMERIC,
    rating_count       INTEGER,
    reference          TEXT,
    source_id          INTEGER,
    CONSTRAINT ck_pains       CHECK (pains IS NULL OR pains IN ('Yes', 'No')),
    CONSTRAINT ck_toxicophore CHECK (toxicophore IS NULL OR toxicophore IN ('Yes', 'No')),
    CONSTRAINT fk_unsuitable_compound
        FOREIGN KEY (inchikey) REFERENCES compound (inchikey),
    CONSTRAINT fk_unsuitable_source
        FOREIGN KEY (source_id) REFERENCES bioactivity_source (source_id)
);

CREATE INDEX idx_unsuitable_source ON unsuitable (source_id);


-- what a source says about a compound that has no column of its own: scores,
-- flags, dates, external ids, the names of its control compounds. one property
-- per row so a source can add one without a migration, and ordinal so a
-- property can hold a list. the key makes a reload idempotent.
CREATE TABLE compound_annotation (
    inchikey      VARCHAR(27) NOT NULL,
    source_db     VARCHAR(255) NOT NULL,
    property      VARCHAR(100) NOT NULL,
    ordinal       INTEGER NOT NULL DEFAULT 0,
    value         TEXT,
    CONSTRAINT pk_compound_annotation
        PRIMARY KEY (inchikey, source_db, property, ordinal),
    CONSTRAINT fk_ca_compound
        FOREIGN KEY (inchikey) REFERENCES compound (inchikey)
);


-- the same, for what a source says about a target that is not its composition.
-- keyed on target_id and not on an accession: one accession carries more than
-- one class when two sources, or two records of one source, disagree.
CREATE TABLE target_annotation (
    target_id     INTEGER NOT NULL,
    source_db     VARCHAR(255) NOT NULL,
    property      VARCHAR(100) NOT NULL,
    ordinal       INTEGER NOT NULL DEFAULT 0,
    value         TEXT,
    CONSTRAINT pk_target_annotation
        PRIMARY KEY (target_id, source_db, property, ordinal),
    CONSTRAINT fk_ta_target
        FOREIGN KEY (target_id) REFERENCES target (target_id)
);


-- a dose given to an animal is not a potency: no endpoint, nothing to compare,
-- and no target, so it cannot be a bioactivity row. one row per dose, and
-- dose_raw keeps the string it was read out of.
CREATE TABLE in_vivo (
    id            SERIAL PRIMARY KEY,
    inchikey      VARCHAR(27) NOT NULL,
    organism      VARCHAR(100),
    dose_value    NUMERIC,
    dose_unit     VARCHAR(50),
    route         VARCHAR(20),
    dose_raw      VARCHAR(255),
    source_id     INTEGER,
    CONSTRAINT ck_route
        CHECK (route IS NULL OR route IN ('IV', 'PO', 'IP', 'SC', 'topical')),
    CONSTRAINT fk_iv_compound
        FOREIGN KEY (inchikey) REFERENCES compound (inchikey),
    CONSTRAINT fk_iv_source
        FOREIGN KEY (source_id) REFERENCES bioactivity_source (source_id)
);


-- the literature a source cites for a compound, where it cites a reading list
-- rather than a paper per measurement. xref_id + source_xref resolves it, the
-- same convention bioactivity_source uses.
CREATE TABLE compound_reference (
    id            SERIAL PRIMARY KEY,
    inchikey      VARCHAR(27) NOT NULL,
    xref_id       VARCHAR(255),
    source_xref   VARCHAR(400),
    raw           TEXT,
    CONSTRAINT fk_cr_compound
        FOREIGN KEY (inchikey) REFERENCES compound (inchikey)
);


-- a number a preprocessor could read but not trust: a factor of 10^n it refuses
-- to apply, a reciprocal rate constant, a range that reads high to low. kept out
-- of bioactivity on purpose, and kept here rather than dropped so a curator can
-- see them. reason says why, raw is the string it came from.
CREATE TABLE quarantine (
    id                SERIAL PRIMARY KEY,
    inchikey          VARCHAR(27) NOT NULL,
    target_id         INTEGER,
    reason            VARCHAR(255) NOT NULL,
    raw               TEXT,
    fragment          TEXT,
    relation          VARCHAR(5),
    value             NUMERIC,
    unit              VARCHAR(50),
    bioactivity_type  VARCHAR(50),
    assay_description TEXT,
    source_id         INTEGER,
    CONSTRAINT fk_q_compound
        FOREIGN KEY (inchikey) REFERENCES compound (inchikey),
    CONSTRAINT fk_q_target
        FOREIGN KEY (target_id) REFERENCES target (target_id),
    CONSTRAINT fk_q_source
        FOREIGN KEY (source_id) REFERENCES bioactivity_source (source_id)
);


-- a record a source published that could not become a compound at all, almost
-- always because it has no structure and the InChIKey is the key everywhere.
-- the one table with no foreign key, because its rows are about absence.
CREATE TABLE skipped_compound (
    id            SERIAL PRIMARY KEY,
    name          VARCHAR(255) NOT NULL,
    source_db     VARCHAR(255) NOT NULL,
    reason        VARCHAR(255),
    portal_path   VARCHAR(255),
    targets       INTEGER,
    validations   INTEGER
);

CREATE INDEX idx_compound_annotation_property ON compound_annotation (property);
CREATE INDEX idx_in_vivo_compound             ON in_vivo (inchikey);
CREATE INDEX idx_compound_reference_compound  ON compound_reference (inchikey);
CREATE INDEX idx_quarantine_compound          ON quarantine (inchikey);
