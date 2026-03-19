#pragma once

// clang-format off
#define GETSET( NAME, CLASS, TYPE, GETFUNCNAME, SETFUNCNAME, GETSTR, SETSTR )                                        \
  .def_prop_rw(                                                                                                      \
      #NAME, []( const CLASS &self ) { return self.GETFUNCNAME(); },                                                 \
      []( CLASS &self, TYPE &val ) { self.SETFUNCNAME( val ); }, nb::for_getter( GETSTR ), nb::for_setter( SETSTR ), \
      nb::for_setter( nb::arg( #NAME ) ) )

#define PDB_WRITE_OPTS(X) \
  X(minimal_file)   \
  X(atom_records)   \
  X(seqres_records)   \
  X(ssbond_records)   \
  X(link_records)   \
  X(cispep_records)   \
  X(cryst1_record)    \
  X(ter_records)    \
  X(conect_records)   \
  X(end_record)   \
  X(numbered_ter)   \
  X(ter_ignores_type)   \
  X(use_linkr)    \
  X(preserve_serial)

#define GROUPS_FIELDS(X)  \
    X(atoms)            \
    X(block_name)       \
    X(entry)            \
    X(database_status)  \
    X(author)           \
    X(cell)             \
    X(symmetry)         \
    X(entity)           \
    X(entity_poly)      \
    X(struct_ref)       \
    X(chem_comp)        \
    X(exptl)            \
    X(diffrn)           \
    X(reflns)           \
    X(refine)           \
    X(title_keywords)   \
    X(ncs)              \
    X(struct_asym)      \
    X(origx)            \
    X(struct_conf)      \
    X(struct_sheet)     \
    X(struct_biol)      \
    X(assembly)         \
    X(conn)             \
    X(cis)              \
    X(modres)           \
    X(scale)            \
    X(atom_type)        \
    X(entity_poly_seq)  \
    X(tls)              \
    X(software)         \
    X(group_pdb)        \
    X(auth_all)         \

#define SETOPTS( OPTS, ITEM, FIELD, MATCHED )          \
  do {                             \
    std::string key = nb::cast<std::string>(ITEM.first); \
    bool val = nb::cast<bool>(ITEM.second); \
    if (key == #FIELD) { opts.FIELD = val; MATCHED = true; } \
  } while (0)
// clang-format on
