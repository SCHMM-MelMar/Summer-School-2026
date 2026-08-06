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
