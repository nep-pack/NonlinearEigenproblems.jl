######################################################################################
# Compares two benchmark JSON files that were previously created by run_benchmark.jl
######################################################################################
include("benchmark_utils.jl")

const WHITE = 15

color(c, text) = @sprintf("\e[38;5;%dm%s\e[0m", c, text)

const TIME_TOLERANCE_PERCENT = [
    (-10, 10)
    ( -5, 114)
    (  5, 7)
    ( 10, 214)
    (Inf, 9)]

time_color(diff) =
    TIME_TOLERANCE_PERCENT[findfirst(t -> (isnan(diff) ? 0 : diff) < t[1], TIME_TOLERANCE_PERCENT)][2]

const MEMORY_TOLERANCE_PERCENT = [
    ( -5, 10)
    ( -1, 114)
    (  1, 7)
    (  5, 214)
    (Inf, 9)]

memory_color(diff) =
    MEMORY_TOLERANCE_PERCENT[findfirst(t -> (isnan(diff) ? 0 : diff) < t[1], MEMORY_TOLERANCE_PERCENT)][2]

function print_table_header(test_hierarchy, time_memory_width, time_padding, memory_padding)
    column_width = recursive_mapreduce(x -> length(x.indented_name), max, parent -> values(parent.children), test_hierarchy)
    column_title = "Test name"
    divider = color(WHITE, "─"^(column_width + time_memory_width)) * "\n"
    header = color(WHITE, @sprintf("%s%sTime%sMemory", column_title,
        " "^(column_width - length(column_title) + time_padding), " "^memory_padding))
    @printf("%s%s\n%s", divider, header, divider)
    return column_width, divider
end

print_benchmark_details(title, d) =
    println("\n",
        color(WHITE, title), "\n",
        color(WHITE, "Run:"), " $(d["time"]); $(d["config"])\n",
        color(WHITE, "CPU:"), " $(d["cpu"])\n",
        color(WHITE, "Git:"), " $(d["git"])")

"Prints a single JSON file benchmark to stdout."
function print_benchmark(baseline_file)
    base_dict = load_benchmark_from_file(baseline_file)
    base = base_dict["benchmark"][1]

    root = get_testset_hierarchy(base)
    column_width, divider = print_table_header(root, 23, 6, 7)

    iterate_hierarchy(root) do test
        trial = test.trial[2]
        @printf("%s%s    %s     %s\n",
            test.indented_name,
            " "^(column_width - length(test.indented_name)),
            prettytime(trial.time, true), prettymemory(trial.memory, true))
    end

    print(divider)

    print_benchmark_details("Details", base_dict)
end

"Prints the comparison of two JSON file benchmarks to stdout."
function print_benchmark_comparison(baseline_file, alternative_file)
    base_dict = load_benchmark_from_file(baseline_file)
    base = base_dict["benchmark"][1]

    alt_dict = load_benchmark_from_file(alternative_file)
    alt = alt_dict["benchmark"][1]

    root = get_testset_hierarchy(alt)
    column_width, divider = print_table_header(root, 39, 9, 15)

    iterate_hierarchy(root) do test
        key = test.trial[1]
        trial = test.trial[2]

        if haskey(base, key)
            j = judge(trial, base[key])
            time_diff = 100 * (j.ratio.time - 1)
            memory_diff = 100 * (j.ratio.memory - 1)
        else
            time_diff = NaN
            memory_diff = NaN
        end

        @printf("%s%s    %s     %s\n",
            test.indented_name,
            " "^(column_width - length(test.indented_name)),
            color(time_color(time_diff), prettytime(trial.time, true) * " (" * prettypercent(time_diff, true) * ")"),
            color(memory_color(memory_diff), prettymemory(trial.memory, true) * " (" * prettypercent(memory_diff, true) * ")"))
    end

    print(divider)

    print_benchmark_details("Baseline", base_dict)
    print_benchmark_details("Alternative", alt_dict)
end

if !(1 <= length(ARGS) <= 2)
    println("Specify one or two benchmark JSON files as input:\n\n$(basename(@__FILE__)) baseline_file [alternative_file]")
else
    length(ARGS) > 1 ? print_benchmark_comparison(ARGS[1], ARGS[2]) : print_benchmark(ARGS[1])
end
