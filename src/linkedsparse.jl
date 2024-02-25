import Base: zero, iszero, iterate

mutable struct LinkedSparseEntry{T}
    const value::BasisVectorEntry{T}
    next::Union{Nothing,LinkedSparseEntry{T}}
end

mutable struct LinkedSparseBasisVector{T} <: AbstractBasisVector{T}
    # this single linked list implementation can be optimized by improving caching and allocation on modern cpus
    # by storing adjacent nodes/entries in continous memory, e.g. by having a vector of nodes with insertions/deletions only at the back
    # and by keeping entries on deletions for reuse rather than freeing them immediately to prevent allocations at the cost of more memory usage  
    # or by preallocating entries for later use
    const basis::Basis{BasisElement{T}}
    first::Union{Nothing,LinkedSparseEntry{T}}

    LinkedSparseBasisVector(basis::Basis{BasisElement{T}}) where {T} = new{T}(basis, nothing)

    LinkedSparseBasisVector(b::BasisElement{T}) where {T} = new{T}(basis(b), LinkedSparseEntry(BasisVectorEntry(b, one(base_ring(basis(b)))), nothing))

    LinkedSparseBasisVector(b::BasisElement{T}, c::T) where {T} = new{T}(basis(b), iszero(c) ? nothing : LinkedSparseEntry(BasisVectorEntry(b, c), nothing))
end

basis(v::LinkedSparseBasisVector) = v.basis

zero(v::LinkedSparseBasisVector) = LinkedSparseBasisVector(basis(v))

function zero!(v::LinkedSparseBasisVector)
    v.first = nothing

    return v
end

function deleteafter!(::LinkedSparseBasisVector{T}, e::LinkedSparseEntry{T}) where {T}
    d = e.next.next
    e.next = d
    return d  # return element now in the position of deleted element, i.e. new element after e
end

function deleteafter!(v::LinkedSparseBasisVector{T}, ::Nothing) where {T}
    d = v.first.next
    v.first = d
    return d  # return element now in the position of deleted element, i.e. new first element
end

function insertafter!(::LinkedSparseBasisVector{T}, e::LinkedSparseEntry{T}, elem::BasisElement{T}, coeff::T) where {T}
    d = LinkedSparseEntry(BasisVectorEntry(elem, coeff), e.next)
    e.next = d
    return d  # return new element
end

function insertafter!(v::LinkedSparseBasisVector{T}, ::Nothing, elem::BasisElement{T}, coeff::T) where {T}
    d = LinkedSparseEntry(BasisVectorEntry(elem, coeff), v.first)
    v.first = d
    return d  # return new element
end

function add_scaled!(v₁::LinkedSparseBasisVector{T}, v₂::LinkedSparseBasisVector{T}, c::T, last₁=nothing) where {T}
    @assert v₁ !== v₂

    iszero(c) && return v₁  # if c = 0 no work needs to be done

    current₁ = last₁ === nothing ? v₁.first : last₁.next  # start at first or at specified position (unsafe)
    current₂ = v₂.first  # start at first
    last₂ = nothing

    while current₁ !== nothing && current₂ !== nothing
        val₂ = current₂.value
        b₂ = val₂.elem

        if isdeleted(b₂)
            add_scaled!(v₂, replacement(b₂), val₂.coeff, current₂)  # add replacement of current₂ (recursively)
            current₂ = deleteafter!(v₂, last₂)  # delete current₂ and go to next on v₂, thus last₂ does not change
        else
            val₁ = current₁.value
            b₁ = val₁.elem

            if b₁ < b₂
                # insert new entry between last₁ and current₁ and stay at position on v₁, thus last₁ changes but current₁ does not
                last₁ = insertafter!(v₁, last₁, b₂, val₂.coeff * c)  # not zero as fields are integral domains

                # go to next on v₂
                last₂ = current₂
                current₂ = current₂.next
            elseif b₁ > b₂
                # go to next on v₁
                last₁ = current₁
                current₁ = current₁.next
            else  # b₁ == b₂
                val₁.coeff += val₂.coeff * c
                if iszero(val₁.coeff)
                    current₁ = deleteafter!(v₁, last₁)  # delete current₁ and go to next on v₁, thus last₁ does not change
                else
                    # go to next on v₁
                    last₁ = current₁
                    current₁ = current₁.next
                end

                # go to next on v₂
                last₂ = current₂
                current₂ = current₂.next
            end
        end
    end

    while current₂ !== nothing
        val₂ = current₂.value
        b₂ = val₂.elem

        if isdeleted(b₂)
            add_scaled!(v₂, replacement(b₂), val₂.coeff, current₂)  # add replacement of current₂ (recursively)
            current₂ = deleteafter!(v₂, last₂)  # delete current₂ and go to next on v₂, thus last₂ does not change
        else
            # to get here current₁ must be nothing and last₁ was the last entry on v₁
            # insert new entry after last₁ and stay at end position on v₁
            last₁ = insertafter!(v₁, last₁, b₂, val₂.coeff * c)  # not zero as fields are integral domains

            # go to next on v₂
            last₂ = current₂
            current₂ = current₂.next
        end
    end

    return v₁
end


function add!(v₁::LinkedSparseBasisVector{T}, v₂::LinkedSparseBasisVector{T}) where {T}
    @assert v₁ !== v₂

    current₁ = v₁.first  # start at first
    current₂ = v₂.first  # start at first
    last₁ = nothing
    last₂ = nothing

    while current₁ !== nothing && current₂ !== nothing
        val₂ = current₂.value
        b₂ = val₂.elem

        if isdeleted(b₂)
            add_scaled!(v₂, replacement(b₂), val₂.coeff, current₂)  # add replacement of current₂
            current₂ = deleteafter!(v₂, last₂)  # delete current₂ and go to next on v₂, thus last₂ does not change
        else
            val₁ = current₁.value
            b₁ = val₁.elem

            if b₁ < b₂
                # insert new entry between last₁ and current₁ and stay at position on v₁, thus last₁ changes but current₁ does not
                last₁ = insertafter!(v₁, last₁, b₂, val₂.coeff)

                # go to next on v₂
                last₂ = current₂
                current₂ = current₂.next
            elseif b₁ > b₂
                # go to next on v₁
                last₁ = current₁
                current₁ = current₁.next
            else  # b₁ == b₂
                val₁.coeff += val₂.coeff
                if iszero(val₁.coeff)
                    current₁ = deleteafter!(v₁, last₁)  # delete current₁ and go to next on v₁, thus last₁ does not change
                else
                    # go to next on v₁
                    last₁ = current₁
                    current₁ = current₁.next
                end

                # go to next on v₂
                last₂ = current₂
                current₂ = current₂.next
            end
        end
    end

    while current₂ !== nothing
        val₂ = current₂.value
        b₂ = val₂.elem

        if isdeleted(b₂)
            add_scaled!(v₂, replacement(b₂), val₂.coeff, current₂)  # add replacement of current₂
            current₂ = deleteafter!(v₂, last₂)  # delete current₂ and go to next on v₂, thus last₂ does not change
        else
            # to get here current₁ must be nothing and last₁ was the last entry on v₁
            # insert new entry after last₁ and stay at end position on v₁
            last₁ = insertafter!(v₁, last₁, b₂, val₂.coeff)

            # go to next on v₂
            last₂ = current₂
            current₂ = current₂.next
        end
    end

    return v₁
end


function subtract!(v₁::LinkedSparseBasisVector{T}, v₂::LinkedSparseBasisVector{T}) where {T}
    @assert v₁ !== v₂

    current₁ = v₁.first  # start at first
    current₂ = v₂.first  # start at first
    last₁ = nothing
    last₂ = nothing

    while current₁ !== nothing && current₂ !== nothing
        val₂ = current₂.value
        b₂ = val₂.elem

        if isdeleted(b₂)
            add_scaled!(v₂, replacement(b₂), val₂.coeff, current₂)  # add replacement of current₂
            current₂ = deleteafter!(v₂, last₂)  # delete current₂ and go to next on v₂, thus last₂ does not change
        else
            val₁ = current₁.value
            b₁ = val₁.elem

            if b₁ < b₂
                # insert new entry between last₁ and current₁ and stay at position on v₁, thus last₁ changes but current₁ does not
                last₁ = insertafter!(v₁, last₁, b₂, -val₂.coeff)

                # go to next on v₂
                last₂ = current₂
                current₂ = current₂.next
            elseif b₁ > b₂
                # go to next on v₁
                last₁ = current₁
                current₁ = current₁.next
            else  # b₁ == b₂
                val₁.coeff -= val₂.coeff
                if iszero(val₁.coeff)
                    current₁ = deleteafter!(v₁, last₁)  # delete current₁ and go to next on v₁, thus last₁ does not change
                else
                    # go to next on v₁
                    last₁ = current₁
                    current₁ = current₁.next
                end

                # go to next on v₂
                last₂ = current₂
                current₂ = current₂.next
            end
        end
    end

    while current₂ !== nothing
        val₂ = current₂.value
        b₂ = val₂.elem

        if isdeleted(b₂)
            add_scaled!(v₂, replacement(b₂), val₂.coeff, current₂)  # add replacement of current₂
            current₂ = deleteafter!(v₂, last₂)  # delete current₂ and go to next on v₂, thus last₂ does not change
        else
            # to get here current₁ must be nothing and last₁ was the last entry on v₁
            # insert new entry after last₁ and stay at end position on v₁
            last₁ = insertafter!(v₁, last₁, b₂, -val₂.coeff)

            # go to next on v₂
            last₂ = current₂
            current₂ = current₂.next
        end
    end

    return v₁
end

# find the next undeleted entry of v after last
function next(v::LinkedSparseBasisVector{T}, last::Union{Nothing,LinkedSparseEntry{T}}) where {T}
    current = last === nothing ? v.first : last.next  # start at first or where we left off

    while current !== nothing
        val = current.value
        b = val.elem

        isdeleted(b) || return current  # found undeleted entry
        add_scaled!(v, replacement(b), val.coeff, current)  # add replacement of current
        current = deleteafter!(v, last)  # delete current and go to next, thus last does not change
    end

    return nothing  # found no undeleted entry
end

# the iterate functions must be short, so they are inlined by the compiler and no tuples for the return values are allocated
iterate(v::LinkedSparseBasisVector) = begin
    e = next(v, nothing)
    e === nothing ? nothing : (e.value, e)
end

iterate(v::LinkedSparseBasisVector{T}, last::LinkedSparseEntry{T}) where {T} = begin
    e = next(v, last)
    e === nothing ? nothing : (e.value, e)
end

last(v::LinkedSparseBasisVector) = begin
    e = next(v, nothing)
    e === nothing ? nothing : e.value
end

iszero(v::LinkedSparseBasisVector) = next(v, nothing) === nothing


function scale!(v::LinkedSparseBasisVector{T}, c::T) where {T}
    iszero(c) && return zero!(v)  # if c = 0 just set v = 0

    current = v.first  # start at first
    last = nothing

    while current !== nothing
        val = current.value
        b = val.elem

        if isdeleted(b)
            add_scaled!(v, replacement(b), val.coeff, current)  # add replacement of current
            current = deleteafter!(v, last)  # delete current and go to next, thus last does not change
        else
            val.coeff *= c  # not zero as fields are integral domains

            # go to next
            last = current
            current = current.next
        end
    end

    return v
end
