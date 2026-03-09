module RavenTemplateTODOs

export RAVEN_TEMPLATE_MODEL_STATUS,
       build_raven_ubcwm,
       build_raven_hbv_ec,
       build_raven_canadian_shield,
       build_raven_mohyse,
       build_raven_hmets,
       build_raven_hypr,
       build_raven_awbm,
       build_raven_sac_sma,
       build_raven_routing_only,
       build_raven_blended_v1,
       build_raven_blended_v2

const RAVEN_TEMPLATE_MODEL_STATUS = (
    ubcwm = (status = :todo, note = "TODO: requires fuller UBCWM precipitation-energy and routing support."),
    hbv_ec = (status = :todo, note = "TODO: requires HBV-EC process sequence verification from the manual."),
    hbv_light = (status = :implemented, note = "Implemented in `hbv_light.jl`."),
    gr4j = (status = :implemented, note = "Implemented in `gr4j.jl`."),
    canadian_shield = (status = :todo, note = "TODO: requires Appendix F Canadian Shield process ordering and wetland routing assembly."),
    mohyse = (status = :todo, note = "TODO: missing exact MOHYSE routing and store transfer sequence."),
    hmets = (status = :todo, note = "TODO: missing exact HMETS partition and routing assembly."),
    hypr = (status = :todo, note = "TODO: missing HYPR routing and store-coupling details."),
    hymod = (status = :implemented, note = "Implemented in `hymod.jl`."),
    awbm = (status = :todo, note = "TODO: missing AWBM area-fraction store family."),
    sac_sma = (status = :todo, note = "TODO: missing full SAC-SMA tension/free-water process family."),
    routing_only = (status = :todo, note = "TODO: needs routing-only wrapper around routing components."),
    blended_v1 = (status = :todo, note = "TODO: requires explicit blended-model process ordering and calibration parameter mapping."),
    blended_v2 = (status = :todo, note = "TODO: requires explicit blended-model process ordering and calibration parameter mapping."),
)

_todo(name::AbstractString, note::AbstractString) =
    throw(ArgumentError("$(name) is not implemented yet. $(note)"))

build_raven_ubcwm() = _todo("raven_ubcwm", RAVEN_TEMPLATE_MODEL_STATUS.ubcwm.note)
build_raven_hbv_ec() = _todo("raven_hbv_ec", RAVEN_TEMPLATE_MODEL_STATUS.hbv_ec.note)
build_raven_canadian_shield() = _todo("raven_canadian_shield", RAVEN_TEMPLATE_MODEL_STATUS.canadian_shield.note)
build_raven_mohyse() = _todo("raven_mohyse", RAVEN_TEMPLATE_MODEL_STATUS.mohyse.note)
build_raven_hmets() = _todo("raven_hmets", RAVEN_TEMPLATE_MODEL_STATUS.hmets.note)
build_raven_hypr() = _todo("raven_hypr", RAVEN_TEMPLATE_MODEL_STATUS.hypr.note)
build_raven_awbm() = _todo("raven_awbm", RAVEN_TEMPLATE_MODEL_STATUS.awbm.note)
build_raven_sac_sma() = _todo("raven_sac_sma", RAVEN_TEMPLATE_MODEL_STATUS.sac_sma.note)
build_raven_routing_only() = _todo("raven_routing_only", RAVEN_TEMPLATE_MODEL_STATUS.routing_only.note)
build_raven_blended_v1() = _todo("raven_blended_v1", RAVEN_TEMPLATE_MODEL_STATUS.blended_v1.note)
build_raven_blended_v2() = _todo("raven_blended_v2", RAVEN_TEMPLATE_MODEL_STATUS.blended_v2.note)

end

