"Calculate the Kolmogorov length scale `η` from viscousity and average TKE dissapation."
η(ν, ϵ) = (ν^3 / ϵ)^(1/4)

"Return the mean from a `FieldTimeSeries` that is `OnDisk()`."
function field_ts_timemean(field_ts::FieldTimeSeries)

    t = field_ts.times
    field_data = field_ts[1].data
    for i ∈ 2:length(t)
        field_data .+= field_ts[i].data
    end

    return field_data ./ length(t)

end
"Find the minimum `η` from the output `ϵ` time series."
function minimum_η(ϵ::FieldTimeSeries; ν = 1e-6)

    t = ϵ.times
    minimum_η_t = similar(t)
    for i ∈ eachindex(t)
        minimum_η_t[i] = minimum(η.(ν, ϵ[i]))
    end

    return minimum(minimum_η_t)

end
