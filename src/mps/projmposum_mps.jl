mutable struct ProjMPOSum_MPS
    pmo::Vector{ProjMPO}
    pm::Vector{ProjMPS}
    weight::Float64
  end
  
  copy(P::ProjMPOSum_MPS) = ProjMPO_MPS(copy(P.pmo), copy.(P.pm), P.weight)
  
  function ProjMPOSum_MPS(mpov::Vector{MPO}, mpsv::Vector{MPS}; weight=1.0)
    return ProjMPOSum_MPS([ProjMPO(M) for M in mpov], [ProjMPS(m) for m in mpsv], weight)
  end
  
  #Maybe having this make no sense anymore?
  ProjMPO_MPS(mpov::Vector{MPO}, Ms::MPS...; weight=1.0) = ProjMPO_MPS(mpov, [Ms...], weight)
  
  nsite(P::ProjMPOSum_MPS) = nsite(P.pmo[1])
  
  function set_nsite!(Ps::ProjMPOSum_MPS, nsite)
    for P in Ps.pmo
        set_nsite!(P, nsite)
    end
    for P in Ps.pm
      set_nsite!(P, nsite)
    end
    return Ps
  end
  
  Base.length(P::ProjMPOSum_MPS) = length(P.pmo[1])
  
  function product(P::ProjMPOSum_MPS, v::ITensor)::ITensor
    Pv = product(P.pmo[1], v)
    for n in 2:length(P.pmo)
        Pv += product(P.pmo[n], v)
    end
    for p in P.pm
      Pv += P.weight * product(p, v)
    end
    return Pv
  end
  
  function Base.eltype(P::ProjMPOSum_MPS)
    elT = eltype(P.pmo[1])
    for n in 2:length(P.pmo)
        elT = promote_type(elT, eltype(P.pmo[n]))
    end
    for p in P.pm
      elT = promote_type(elT, eltype(p))
    end
    return elT
  end
  
  (P::ProjMPOSum_MPS)(v::ITensor) = product(P, v)
  
  Base.size(P::ProjMPOSum_MPS) = size(P.pmo[1])
  
  function position!(P::ProjMPOSum_MPS, psi::MPS, pos::Int)
    for M in P.pmo
        position!(M, psi, pos)
    end
    for p in P.pm
      position!(p, psi, pos)
    end
    return P
  end
  

  function noiseterm(P::ProjMPOSum_MPS, phi::ITensor, dir::String)
    nt = noiseterm(P.pmo[1], phi, dir)
    for n in 2:length(P.pmo)
      nt += noiseterm(P.pmo[n], phi, dir)
    end
    return nt
  end

