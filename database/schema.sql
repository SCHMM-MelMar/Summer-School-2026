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


-- a compound set is a collection of compounds, and the loader records one per
-- staging directory, named after it, so every compound says which source it
-- came in from. category is what kind of collection it is, left at 'library'
-- unless somebody curates it. a compound is in as many sets as claim it, which
-- is the point: the overlap between two sources is a join and not a manual
-- comparison, and that is why membership is its own table rather than a column
-- on the compound.
CREATE TABLE compound_set (
    set_id        SERIAL PRIMARY KEY,
    name          VARCHAR(255) NOT NULL,
    category      VARCHAR(50) NOT NULL,
    source_db     VARCHAR(255),
    description   TEXT,
    CONSTRAINT uq_compound_set UNIQUE (name),
    CONSTRAINT ck_set_category
        CHECK (category IN ('chemical probe', 'chemogenomic', 'drug',
                            'negative control', 'library', 'other'))
);


CREATE TABLE compound_set_member (
    set_id        INTEGER NOT NULL,
    inchikey      VARCHAR(27) NOT NULL,
    CONSTRAINT pk_compound_set_member PRIMARY KEY (set_id, inchikey),
    CONSTRAINT fk_csm_set
        FOREIGN KEY (set_id) REFERENCES compound_set (set_id),
    CONSTRAINT fk_csm_compound
        FOREIGN KEY (inchikey) REFERENCES compound (inchikey)
);


-- superfamily is UniProt's own classification of the protein, read from
-- reference/uniprot_protein_families.tsv rather than from any one source, so
-- every source agrees on it. it is the whole string UniProt gives, outermost
-- group first: "Small GTPase superfamily, Ras family". not every entry has one.
-- this is a property of one protein, not of a target: a complex has as many
-- classifications as it has members, and a target of type 'family' is our own
-- curated grouping, which is a different thing entirely.
CREATE TABLE uniprot (
    uniprot_id    VARCHAR(20) PRIMARY KEY,
    entrez_gene   VARCHAR(50),
    hgnc          VARCHAR(50),
    species       VARCHAR(100),
    superfamily   VARCHAR(255)
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
CREATE INDEX idx_set_member_compound   ON compound_set_member (inchikey);
CREATE INDEX idx_bioactivity_compound  ON bioactivity (inchikey);
CREATE INDEX idx_bioactivity_target    ON bioactivity (target_id);
CREATE INDEX idx_bioactivity_source    ON bioactivity (source_id);
CREATE INDEX idx_bioactivity_scope     ON bioactivity (bioactivity_type, unit);
CREATE INDEX idx_compound_name         ON compound (name);
CREATE INDEX idx_chembl_inchikey       ON chembl (inchikey);


-- bioactivity stores source_id, an integer, because the name of the source
-- lives in one place and not on every measurement. this view is the same rows
-- with the compound, the target and the source spelled out, so reading does
-- not mean remembering three joins. db.bioactivities() is this view filtered.
CREATE VIEW bioactivity_flat AS
SELECT b.id, b.inchikey, c.name AS compound, b.target_id,
       t.type AS target_type, t.name AS target, b.moa,
       b.bioactivity_type, b.relation, b.value, b.unit,
       b.assay_type, b.assay_description, b.cell_line,
       b.concentration, b.concentration_unit,
       b.source_id, s.source_db, s.source, b.source_xref,
       CASE WHEN s.xref_id IS NOT NULL AND b.source_xref IS NOT NULL
            THEN s.xref_id || b.source_xref END AS source_url
  FROM bioactivity b
  JOIN compound c ON c.inchikey = b.inchikey
  JOIN target t ON t.target_id = b.target_id
  LEFT JOIN bioactivity_source s ON s.source_id = b.source_id;


-- membership is a link table, so this view is one row per compound and set: a
-- compound in three libraries is three rows, and a compound in none still
-- shows up with the set columns empty. it is the join surface, the thing you
-- merge your own frame onto. db.compounds() is the same data rolled up to one
-- row per compound instead.
CREATE VIEW compound_flat AS
SELECT c.inchikey, c.name, c.smiles,
       s.set_id, s.name AS set_name, s.category, s.source_db
  FROM compound c
  LEFT JOIN compound_set_member m ON m.inchikey = c.inchikey
  LEFT JOIN compound_set s ON s.set_id = m.set_id;


-- target and uniprot have no column in common, the link is target_uniprot.
-- this view does both joins so you get one flat frame instead of three.
-- LEFT so a target with no accession yet still shows up.
CREATE VIEW target_flat AS
SELECT t.target_id, t.type, t.name,
       tu.uniprot_id, u.hgnc, u.species, u.entrez_gene, u.superfamily
  FROM target t
  LEFT JOIN target_uniprot tu ON tu.target_id = t.target_id
  LEFT JOIN uniprot u ON u.uniprot_id = tu.uniprot_id;
