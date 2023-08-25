"Calculate the Kolmogorov length scale `η` from viscousity and average TKE dissapation."
η(ν, ϵ) = (ν^3 / ϵ)^(1/4)

function field_ts_mean(field_ts::FieldTimeSeries)

    t = field_ts.times
    field_data = field_ts[1].data
    for i ∈ 2:length(t)
        field_data .+= field_ts[i].data
    end

    return field_data ./ length(t)

end
