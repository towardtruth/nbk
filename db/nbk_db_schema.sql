--
-- PostgreSQL database dump
--

-- Dumped from database version 12.3
-- Dumped by pg_dump version 12.3

-- Started on 2021-05-18 09:59:08

SET statement_timeout = 0;
SET lock_timeout = 0;
SET idle_in_transaction_session_timeout = 0;
SET client_encoding = 'UTF8';
SET standard_conforming_strings = on;
SELECT pg_catalog.set_config('search_path', '', false);
SET check_function_bodies = false;
SET xmloption = content;
SET client_min_messages = warning;
SET row_security = off;

SET default_tablespace = '';

SET default_table_access_method = heap;

--
-- TOC entry 202 (class 1259 OID 16411)
-- Name: cutoffs; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.cutoffs (
    coid integer NOT NULL,
    gp_type character(1) NOT NULL,
    tsid integer NOT NULL,
    percentile real NOT NULL,
    cutoff double precision NOT NULL
);


--
-- TOC entry 203 (class 1259 OID 16414)
-- Name: cutoffs_coid_seq; Type: SEQUENCE; Schema: public; Owner: -
--

CREATE SEQUENCE public.cutoffs_coid_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


--
-- TOC entry 2922 (class 0 OID 0)
-- Dependencies: 203
-- Name: cutoffs_coid_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: -
--

ALTER SEQUENCE public.cutoffs_coid_seq OWNED BY public.cutoffs.coid;


--
-- TOC entry 204 (class 1259 OID 16424)
-- Name: feature_set; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.feature_set (
    fsid integer NOT NULL,
    fs_name character(32),
    label_data_type character(2),
    description character varying(1024),
    class_size smallint
);


--
-- TOC entry 205 (class 1259 OID 16430)
-- Name: feature_set_fsid_seq; Type: SEQUENCE; Schema: public; Owner: -
--

CREATE SEQUENCE public.feature_set_fsid_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


--
-- TOC entry 2923 (class 0 OID 0)
-- Dependencies: 205
-- Name: feature_set_fsid_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: -
--

ALTER SEQUENCE public.feature_set_fsid_seq OWNED BY public.feature_set.fsid;


--
-- TOC entry 206 (class 1259 OID 16432)
-- Name: features; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.features (
    fid integer NOT NULL,
    fsid integer NOT NULL,
    feature_name character(32),
    description character(1024),
    assigned_genes text,
    corresp_tissue integer,
    feature_size integer
);


--
-- TOC entry 207 (class 1259 OID 16438)
-- Name: features_fid_seq; Type: SEQUENCE; Schema: public; Owner: -
--

CREATE SEQUENCE public.features_fid_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


--
-- TOC entry 2924 (class 0 OID 0)
-- Dependencies: 207
-- Name: features_fid_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: -
--

ALTER SEQUENCE public.features_fid_seq OWNED BY public.features.fid;


--
-- TOC entry 208 (class 1259 OID 16440)
-- Name: gene_ids; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.gene_ids (
    gnid integer NOT NULL,
    gene_id character(20) NOT NULL
);


--
-- TOC entry 209 (class 1259 OID 16443)
-- Name: gene_ids_gnid_seq; Type: SEQUENCE; Schema: public; Owner: -
--

CREATE SEQUENCE public.gene_ids_gnid_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


--
-- TOC entry 2925 (class 0 OID 0)
-- Dependencies: 209
-- Name: gene_ids_gnid_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: -
--

ALTER SEQUENCE public.gene_ids_gnid_seq OWNED BY public.gene_ids.gnid;


--
-- TOC entry 210 (class 1259 OID 16445)
-- Name: gene_sequence_gsid_seq; Type: SEQUENCE; Schema: public; Owner: -
--

CREATE SEQUENCE public.gene_sequence_gsid_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


--
-- TOC entry 211 (class 1259 OID 16447)
-- Name: gene_sequence; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.gene_sequence (
    gsid integer DEFAULT nextval('public.gene_sequence_gsid_seq'::regclass) NOT NULL,
    gnid integer,
    gene_sub_id character(4),
    seq text,
    seq_len integer,
    seq_type character(2)
);


--
-- TOC entry 212 (class 1259 OID 16454)
-- Name: kmers; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.kmers (
    kmid integer NOT NULL,
    k integer,
    gsid integer,
    kmer character varying(10),
    freq integer
);


--
-- TOC entry 213 (class 1259 OID 16457)
-- Name: kmers_kmid_seq; Type: SEQUENCE; Schema: public; Owner: -
--

CREATE SEQUENCE public.kmers_kmid_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


--
-- TOC entry 2926 (class 0 OID 0)
-- Dependencies: 213
-- Name: kmers_kmid_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: -
--

ALTER SEQUENCE public.kmers_kmid_seq OWNED BY public.kmers.kmid;


--
-- TOC entry 214 (class 1259 OID 16459)
-- Name: kmers_promoter; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.kmers_promoter (
    kmid integer NOT NULL,
    gsid integer,
    k integer,
    kmer character varying(20),
    freq integer
);


--
-- TOC entry 215 (class 1259 OID 16462)
-- Name: kmers_promoter_kmid_seq; Type: SEQUENCE; Schema: public; Owner: -
--

CREATE SEQUENCE public.kmers_promoter_kmid_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


--
-- TOC entry 2927 (class 0 OID 0)
-- Dependencies: 215
-- Name: kmers_promoter_kmid_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: -
--

ALTER SEQUENCE public.kmers_promoter_kmid_seq OWNED BY public.kmers_promoter.kmid;


--
-- TOC entry 216 (class 1259 OID 16472)
-- Name: sample_code; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.sample_code (
    scid integer NOT NULL,
    name character varying(30) NOT NULL
);


--
-- TOC entry 217 (class 1259 OID 16475)
-- Name: tissue_code; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.tissue_code (
    tsid integer NOT NULL,
    name character varying(30) NOT NULL
);


--
-- TOC entry 218 (class 1259 OID 16478)
-- Name: tissue_code_tsid_seq; Type: SEQUENCE; Schema: public; Owner: -
--

CREATE SEQUENCE public.tissue_code_tsid_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


--
-- TOC entry 2928 (class 0 OID 0)
-- Dependencies: 218
-- Name: tissue_code_tsid_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: -
--

ALTER SEQUENCE public.tissue_code_tsid_seq OWNED BY public.tissue_code.tsid;


--
-- TOC entry 219 (class 1259 OID 16480)
-- Name: validations_vdid_seq; Type: SEQUENCE; Schema: public; Owner: -
--

CREATE SEQUENCE public.validations_vdid_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    MAXVALUE 2147483647
    CACHE 1;


--
-- TOC entry 220 (class 1259 OID 16489)
-- Name: walley_data; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.walley_data (
    wdid integer NOT NULL,
    scid integer NOT NULL,
    gnid integer,
    value_type character(1),
    ea_value double precision
);


--
-- TOC entry 221 (class 1259 OID 16492)
-- Name: walley_data_wdid_seq; Type: SEQUENCE; Schema: public; Owner: -
--

CREATE SEQUENCE public.walley_data_wdid_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


--
-- TOC entry 2929 (class 0 OID 0)
-- Dependencies: 221
-- Name: walley_data_wdid_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: -
--

ALTER SEQUENCE public.walley_data_wdid_seq OWNED BY public.walley_data.wdid;


--
-- TOC entry 2744 (class 2604 OID 16496)
-- Name: cutoffs coid; Type: DEFAULT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.cutoffs ALTER COLUMN coid SET DEFAULT nextval('public.cutoffs_coid_seq'::regclass);


--
-- TOC entry 2745 (class 2604 OID 16498)
-- Name: feature_set fsid; Type: DEFAULT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.feature_set ALTER COLUMN fsid SET DEFAULT nextval('public.feature_set_fsid_seq'::regclass);


--
-- TOC entry 2746 (class 2604 OID 16499)
-- Name: features fid; Type: DEFAULT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.features ALTER COLUMN fid SET DEFAULT nextval('public.features_fid_seq'::regclass);


--
-- TOC entry 2747 (class 2604 OID 16500)
-- Name: gene_ids gnid; Type: DEFAULT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.gene_ids ALTER COLUMN gnid SET DEFAULT nextval('public.gene_ids_gnid_seq'::regclass);


--
-- TOC entry 2749 (class 2604 OID 16501)
-- Name: kmers kmid; Type: DEFAULT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.kmers ALTER COLUMN kmid SET DEFAULT nextval('public.kmers_kmid_seq'::regclass);


--
-- TOC entry 2750 (class 2604 OID 16502)
-- Name: kmers_promoter kmid; Type: DEFAULT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.kmers_promoter ALTER COLUMN kmid SET DEFAULT nextval('public.kmers_promoter_kmid_seq'::regclass);


--
-- TOC entry 2751 (class 2604 OID 16504)
-- Name: tissue_code tsid; Type: DEFAULT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.tissue_code ALTER COLUMN tsid SET DEFAULT nextval('public.tissue_code_tsid_seq'::regclass);


--
-- TOC entry 2752 (class 2604 OID 16505)
-- Name: walley_data wdid; Type: DEFAULT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.walley_data ALTER COLUMN wdid SET DEFAULT nextval('public.walley_data_wdid_seq'::regclass);


--
-- TOC entry 2754 (class 2606 OID 32831)
-- Name: cutoffs cutoffs_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.cutoffs
    ADD CONSTRAINT cutoffs_pkey PRIMARY KEY (coid);


--
-- TOC entry 2757 (class 2606 OID 32837)
-- Name: feature_set feature_set_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.feature_set
    ADD CONSTRAINT feature_set_pkey PRIMARY KEY (fsid);


--
-- TOC entry 2759 (class 2606 OID 32839)
-- Name: features features_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.features
    ADD CONSTRAINT features_pkey PRIMARY KEY (fid);


--
-- TOC entry 2761 (class 2606 OID 32841)
-- Name: gene_ids gene_ids_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.gene_ids
    ADD CONSTRAINT gene_ids_pkey PRIMARY KEY (gnid);


--
-- TOC entry 2766 (class 2606 OID 32843)
-- Name: gene_sequence gene_sequence_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.gene_sequence
    ADD CONSTRAINT gene_sequence_pkey PRIMARY KEY (gsid);


--
-- TOC entry 2771 (class 2606 OID 32845)
-- Name: kmers kmers_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.kmers
    ADD CONSTRAINT kmers_pkey PRIMARY KEY (kmid);


--
-- TOC entry 2774 (class 2606 OID 32847)
-- Name: kmers_promoter kmers_promoter_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.kmers_promoter
    ADD CONSTRAINT kmers_promoter_pkey PRIMARY KEY (kmid);


--
-- TOC entry 2776 (class 2606 OID 32854)
-- Name: sample_code sample_code_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.sample_code
    ADD CONSTRAINT sample_code_pkey PRIMARY KEY (scid);


--
-- TOC entry 2778 (class 2606 OID 32856)
-- Name: tissue_code tissue_code_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.tissue_code
    ADD CONSTRAINT tissue_code_pkey PRIMARY KEY (tsid);


--
-- TOC entry 2764 (class 2606 OID 32860)
-- Name: gene_ids uq_walley_data_gene_id; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.gene_ids
    ADD CONSTRAINT uq_walley_data_gene_id UNIQUE (gene_id);


--
-- TOC entry 2782 (class 2606 OID 32864)
-- Name: walley_data walley_data_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.walley_data
    ADD CONSTRAINT walley_data_pkey PRIMARY KEY (wdid);


--
-- TOC entry 2762 (class 1259 OID 32869)
-- Name: idx_gene_ids_gene_id; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX idx_gene_ids_gene_id ON public.gene_ids USING btree (gene_id);


--
-- TOC entry 2767 (class 1259 OID 32870)
-- Name: idx_gene_sequence_seq_type_gnid_seq_len; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX idx_gene_sequence_seq_type_gnid_seq_len ON public.gene_sequence USING btree (seq_type, gnid, seq_len);


--
-- TOC entry 2768 (class 1259 OID 32871)
-- Name: idx_kmers_k_gsid_kmer; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX idx_kmers_k_gsid_kmer ON public.kmers USING btree (k, gsid, kmer);


--
-- TOC entry 2769 (class 1259 OID 32872)
-- Name: idx_kmers_k_kmer_gsid; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX idx_kmers_k_kmer_gsid ON public.kmers USING btree (k, kmer, gsid);


--
-- TOC entry 2772 (class 1259 OID 32873)
-- Name: idx_kmers_promoter_gsid_k_kmer; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX idx_kmers_promoter_gsid_k_kmer ON public.kmers_promoter USING btree (gsid, k, kmer);


--
-- TOC entry 2755 (class 1259 OID 32874)
-- Name: idx_pg_type_tsid_percentile; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX idx_pg_type_tsid_percentile ON public.cutoffs USING btree (gp_type, tsid, percentile);


--
-- TOC entry 2779 (class 1259 OID 32876)
-- Name: idx_walley_data_value_type_scid_gnid; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX idx_walley_data_value_type_scid_gnid ON public.walley_data USING btree (value_type, scid, gnid);


--
-- TOC entry 2780 (class 1259 OID 32877)
-- Name: idx_walley_data_vaule_type_gnid_scid; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX idx_walley_data_vaule_type_gnid_scid ON public.walley_data USING btree (value_type, gnid, scid);


--
-- TOC entry 2783 (class 2606 OID 32878)
-- Name: cutoffs cutoffs_tsid_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.cutoffs
    ADD CONSTRAINT cutoffs_tsid_fkey FOREIGN KEY (tsid) REFERENCES public.tissue_code(tsid);


--
-- TOC entry 2784 (class 2606 OID 32888)
-- Name: features features_corresp_tissue_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.features
    ADD CONSTRAINT features_corresp_tissue_fkey FOREIGN KEY (corresp_tissue) REFERENCES public.tissue_code(tsid);


--
-- TOC entry 2785 (class 2606 OID 32893)
-- Name: features features_fsid_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.features
    ADD CONSTRAINT features_fsid_fkey FOREIGN KEY (fsid) REFERENCES public.feature_set(fsid);


--
-- TOC entry 2786 (class 2606 OID 32913)
-- Name: gene_sequence gene_sequence_gnid_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.gene_sequence
    ADD CONSTRAINT gene_sequence_gnid_fkey FOREIGN KEY (gnid) REFERENCES public.gene_ids(gnid);


--
-- TOC entry 2787 (class 2606 OID 32918)
-- Name: kmers kmers_gsid_fKey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.kmers
    ADD CONSTRAINT "kmers_gsid_fKey" FOREIGN KEY (gsid) REFERENCES public.gene_sequence(gsid);


--
-- TOC entry 2788 (class 2606 OID 32923)
-- Name: kmers_promoter kmers_promoter_gsid_fKey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.kmers_promoter
    ADD CONSTRAINT "kmers_promoter_gsid_fKey" FOREIGN KEY (gsid) REFERENCES public.gene_sequence(gsid);


--
-- TOC entry 2789 (class 2606 OID 32928)
-- Name: walley_data scid; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.walley_data
    ADD CONSTRAINT scid FOREIGN KEY (scid) REFERENCES public.sample_code(scid);


--
-- TOC entry 2790 (class 2606 OID 32933)
-- Name: walley_data walley_data_gnid_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.walley_data
    ADD CONSTRAINT walley_data_gnid_fkey FOREIGN KEY (gnid) REFERENCES public.gene_ids(gnid);


-- Completed on 2021-05-18 09:59:08

--
-- PostgreSQL database dump complete
--

