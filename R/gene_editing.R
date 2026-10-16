extract_read_region <- function(
    read_seq,
    cigar,
    ref_start,
    x,
    y,
    include_insertions = TRUE
) {
  read_seq <- toupper(as.character(read_seq))
  ref_start <- as.integer(ref_start)
  x <- as.integer(x)
  y <- as.integer(y)

  # 拆解 CIGAR，例如 3M2I2M2D3M
  ops <- regmatches(
    cigar,
    gregexpr("[0-9]+[MIDNSHP=X]", cigar)
  )[[1]]

  if (length(ops) == 0L || identical(ops, character(0))) {
    stop("无法解析 CIGAR：", cigar)
  }

  op_len <- as.integer(sub("[A-Z=]$", "", ops))
  op_code <- sub("^[0-9]+", "", ops)

  # 检查 read 的比对范围是否覆盖完整的 [x, y]
  ref_consumed <- sum(
    op_len[op_code %in% c("M", "D", "N", "=", "X")]
  )
  ref_end <- ref_start + ref_consumed - 1L

  if (x < ref_start || y > ref_end) {
    return(NA_character_)
  }

  # r：当前参考坐标；q：当前 read 序列坐标
  r <- ref_start
  q <- 1L
  out <- character()

  for (i in seq_along(ops)) {
    n <- op_len[i]
    code <- op_code[i]

    if (code %in% c("M", "=", "X")) {
      # 逐碱基检查是否落在目标参考区间
      for (offset in seq_len(n)) {
        ref_pos <- r + offset - 1L
        read_pos <- q + offset - 1L

        if (ref_pos >= x && ref_pos <= y) {
          out <- c(out, substr(read_seq, read_pos, read_pos))
        }
      }

      r <- r + n
      q <- q + n

    } else if (code == "I") {
      # 插入位于 r-1 与 r 之间。
      # 纳入区间起点前的边界（x-1/x），以及区间内部的插入；
      # 不纳入 y 之后的插入。
      position_after <- r - 1L

      if (include_insertions &&
          position_after >= x - 1L &&
          position_after < y) {
        inserted_seq <- substr(read_seq, q, q + n - 1L)
        out <- c(out, strsplit(inserted_seq, "", fixed = TRUE)[[1]])
      }

      q <- q + n

    } else if (code == "D" || code == "N") {
      # D/N 消耗参考坐标，但不消耗 read 序列。
      # 落在目标区间的每个参考位置用一个横线表示。
      overlap_start <- max(r, x)
      overlap_end <- min(r + n - 1L, y)

      if (overlap_start <= overlap_end) {
        out <- c(out, rep("-", overlap_end - overlap_start + 1L))
      }

      r <- r + n

    } else if (code == "S") {
      # 软剪切碱基在 read 序列中，但没有比对到参考序列
      q <- q + n

    } else if (code %in% c("H", "P")) {
      # H 和 P 不消耗 read 或参考坐标
      next
    }
  }

  paste0(out, collapse = "")
}

extract_aligned_region <- function(
    reference,
    read_seq,
    cigar,
    ref_start,
    x,
    y,
    include_insertions = TRUE
) {
  reference <- toupper(as.character(reference))
  read_seq <- toupper(as.character(read_seq))
  ref_start <- as.integer(ref_start)
  x <- as.integer(x)
  y <- as.integer(y)

  if (x < 1L || y < x || y > nchar(reference)) {
    stop("x 和 y 必须是参考序列上的有效 1-based 坐标，且 x <= y。")
  }

  # 解析 CIGAR，例如 3M2I2M2D3M
  ops <- regmatches(
    cigar,
    gregexpr("[0-9]+[MIDNSHP=X]", cigar)
  )[[1]]

  if (length(ops) == 0L || identical(ops, character(0))) {
    stop("无法解析 CIGAR：", cigar)
  }

  op_len <- as.integer(sub("[A-Z=]$", "", ops))
  op_code <- sub("^[0-9]+", "", ops)

  # 计算比对在参考序列上消耗的长度
  ref_consumed <- sum(op_len[op_code %in% c("M", "D", "N", "=", "X")])
  ref_end <- ref_start + ref_consumed - 1L

  # 当前 read 未覆盖整个请求区间
  if (x < ref_start || y > ref_end) {
    return(list(
      reference_aligned = NA_character_,
      read_aligned = NA_character_
    ))
  }

  r <- ref_start  # 当前参考坐标，1-based
  q <- 1L         # 当前 read 坐标，1-based

  ref_out <- character()
  read_out <- character()

  for (i in seq_along(ops)) {
    n <- op_len[i]
    code <- op_code[i]

    if (code %in% c("M", "=", "X")) {
      # 对齐碱基：两边都取序列
      for (offset in seq_len(n)) {
        ref_pos <- r + offset - 1L
        read_pos <- q + offset - 1L

        if (ref_pos >= x && ref_pos <= y) {
          ref_out <- c(
            ref_out,
            substr(reference, ref_pos, ref_pos)
          )
          read_out <- c(
            read_out,
            substr(read_seq, read_pos, read_pos)
          )
        }
      }

      r <- r + n
      q <- q + n

    } else if (code == "I") {
      # 插入：reference 补横线，read 保留插入序列。
      # position_after 表示插入发生在参考序列哪个碱基之后。
      position_after <- r - 1L

      if (include_insertions &&
          position_after >= x - 1L &&
          position_after < y) {
        inserted_seq <- strsplit(
          substr(read_seq, q, q + n - 1L),
          "",
          fixed = TRUE
        )[[1]]

        ref_out <- c(ref_out, rep("-", n))
        read_out <- c(read_out, inserted_seq)
      }

      q <- q + n

    } else if (code %in% c("D", "N")) {
      # 缺失/跳过：reference 保留碱基，read 补横线
      for (offset in seq_len(n)) {
        ref_pos <- r + offset - 1L

        if (ref_pos >= x && ref_pos <= y) {
          ref_out <- c(
            ref_out,
            substr(reference, ref_pos, ref_pos)
          )
          read_out <- c(read_out, "-")
        }
      }

      r <- r + n

    } else if (code == "S") {
      # 软剪切只消耗 read 序列，不对应参考序列
      q <- q + n

    } else if (code %in% c("H", "P")) {
      # 硬剪切和 padding 不消耗 read 或参考序列
      next
    }
  }

  list(
    reference_aligned = paste0(ref_out, collapse = ""),
    read_aligned = paste0(read_out, collapse = "")
  )
}

process_one_read2 <- function(
    read_id,
    read_seq,
    cigar,
    ref_start,
    reference,
    region_start = 1L,
    region_end = nchar(as.character(reference)),
    min_anchor = 0L
) {
  # reference 可传普通字符字符串，也可传 DNAString
  reference <- as.character(reference)

  ref_start <- as.integer(ref_start)
  region_start <- as.integer(region_start)
  region_end <- as.integer(region_end)
  min_anchor <- as.integer(min_anchor)

  # trim reference and reads
  # reference_trim = unname(substr(reference, region_start, region_end))
  # read_trim = extract_read_region(read_seq = read_seq, cigar = cigar, ref_start = ref_start, x = region_start, y = region_end)
  ref_read = extract_aligned_region(reference = reference, read_seq = read_seq,
                                    cigar = cigar, ref_start = ref_start,
                                    x = region_start, y = region_end,
                                    include_insertions = TRUE)

  # 拆分 CIGAR，例如 10S50M2I40M
  ops <- regmatches(
    cigar,
    gregexpr("[0-9]+[MIDNSHP=X]", cigar)
  )[[1]]

  if (length(ops) == 0L || identical(ops, character(0))) {
    stop("无法解析 CIGAR：", cigar)
  }

  op_len <- as.integer(sub("[A-Z=]$", "", ops))
  op_code <- sub("^[0-9]+", "", ops)

  # M、D、N、=、X 消耗参考序列
  ref_consumed <- sum(op_len[op_code %in% c("M", "D", "N", "=", "X")])

  # ref_consumed = cigarillo::cigar_extent_along_ref(cigar)

  # read 在参考序列上的覆盖范围，坐标为 1-based
  aligned_start <- ref_start
  aligned_end <- ref_start + ref_consumed - 1L

  # 检查是否覆盖目标区间及两侧锚定区
  required_start <- max(1L, region_start - min_anchor)
  required_end <- region_end + min_anchor

  if (aligned_start > required_start || aligned_end < required_end) {
    return(list(
      status = "NO_CALL",
      read_id = read_id,
      combination = "NO_CALL",
      events = data.frame(
        read_id = read_id,
        type = "NO_CALL",
        start = NA_integer_,
        end = NA_integer_,
        position_after = NA_integer_,
        ref = NA_character_,
        alt = NA_character_,
        label = "NO_CALL",

        trim_region_start = region_start,
        trim_region_end = region_end,
        reference_aligned_trim = ref_read$reference_aligned,
        read_aligned_trim = ref_read$read_aligned,

        stringsAsFactors = FALSE
      )
    ))
  }

  # 记录替换碱基的位置和等位序列
  mismatch_pos <- integer()
  mismatch_ref <- character()
  mismatch_alt <- character()

  # 记录插入和缺失事件
  indel_events <- list()
  n_indel <- 0L

  # 当前参考位置和 read/query 位置均用 1-based
  r <- ref_start
  q <- 1L

  for (i in seq_along(ops)) {
    n <- op_len[i]
    code <- op_code[i]

    if (code %in% c("M", "=", "X")) {
      # 逐碱基比较 read 与参考序列
      for (offset in seq_len(n)) {
        current_ref_pos <- r + offset - 1L
        current_read_pos <- q + offset - 1L

        if (current_ref_pos >= region_start &&
            current_ref_pos <= region_end) {
          rb <- substr(reference, current_ref_pos, current_ref_pos)
          qb <- substr(read_seq, current_read_pos, current_read_pos)

          if (rb != qb) {
            # if (!identical(rb, qb)) {
            mismatch_pos <- c(mismatch_pos, current_ref_pos)
            mismatch_ref <- c(mismatch_ref, rb)
            mismatch_alt <- c(mismatch_alt, qb)
          }
        }
      }

      r <- r + n
      q <- q + n

    } else if (code == "I") {
      # 插入发生在当前参考位置之前，即插入在 r-1 与 r 之间。
      # position_after 表示插入发生在参考序列第几个碱基之后。
      position_after <- r - 1L

      if (position_after >= region_start - 1L &&
          position_after <= region_end) {
        inserted_seq <- substr(read_seq, q, q + n - 1L)

        n_indel <- n_indel + 1L
        indel_events[[n_indel]] <- data.frame(
          read_id = read_id, ###################################
          type = "IN",
          start = NA_integer_,
          end = NA_integer_,
          position_after = position_after,
          ref = "-",
          alt = inserted_seq,
          label = paste0("IN_after_", position_after, ":", inserted_seq),

          trim_region_start = region_start,
          trim_region_end = region_end,
          reference_aligned_trim = ref_read$reference_aligned,
          read_aligned_trim = ref_read$read_aligned,

          stringsAsFactors = FALSE
        )
      }

      q <- q + n

    } else if (code == "D") {
      # 删除的参考坐标为 r 到 r+n-1
      deletion_start <- r
      deletion_end <- r + n - 1L

      # 只要缺失与分析区间有重叠，就记录整个连续缺失事件
      if (deletion_end >= region_start &&
          deletion_start <= region_end) {
        deleted_seq <- substr(reference, deletion_start, deletion_end)

        n_indel <- n_indel + 1L
        indel_events[[n_indel]] <- data.frame(
          read_id = read_id, ###################################
          type = "DEL",
          start = deletion_start,
          end = deletion_end,
          position_after = NA_integer_,
          ref = deleted_seq,
          alt = "-",
          label = paste0(
            "DEL_", deletion_start, "_", deletion_end,
            ":", deleted_seq
          ),

          trim_region_start = region_start,
          trim_region_end = region_end,
          reference_aligned_trim = ref_read$reference_aligned,
          read_aligned_trim = ref_read$read_aligned,

          stringsAsFactors = FALSE
        )
      }

      r <- r + n

    } else if (code == "N") {
      # N 表示跳过参考序列；若与目标区间重叠，则无法可靠判读
      skipped_start <- r
      skipped_end <- r + n - 1L

      if (skipped_end >= region_start &&
          skipped_start <= region_end) {
        return(list(
          status = "NO_CALL",
          read_id = read_id,
          combination = "NO_CALL",
          events = data.frame(
            read_id = read_id,
            type = "NO_CALL",
            start = NA_integer_,
            end = NA_integer_,
            position_after = NA_integer_,
            ref = NA_character_,
            alt = NA_character_,
            label = "NO_CALL_N_skip",

            trim_region_start = region_start,
            trim_region_end = region_end,
            reference_aligned_trim = ref_read$reference_aligned,
            read_aligned_trim = ref_read$read_aligned,

            stringsAsFactors = FALSE
          )
        ))
      }

      r <- r + n

    } else if (code %in% c("S", "I")) {
      # S 消耗 read 序列，不消耗参考序列
      q <- q + n

    } else if (code %in% c("H", "P")) {
      # H 和 P 不消耗 read 或参考序列
      next
    }
  }

  # 将相邻的替换碱基合并为 SNP 或 MNP 事件
  substitution_events <- list()

  if (length(mismatch_pos) > 0L) {
    groups <- cumsum(c(TRUE, diff(mismatch_pos) > 1L))
    group_ids <- unique(groups)

    for (g in group_ids) {
      ix <- which(groups == g)

      event_start <- mismatch_pos[ix[1]]
      event_end <- mismatch_pos[ix[length(ix)]]
      ref_allele <- paste0(mismatch_ref[ix], collapse = "")
      alt_allele <- paste0(mismatch_alt[ix], collapse = "")

      event_type <- if (length(ix) == 1L) "SNP" else "MNP"

      substitution_events[[length(substitution_events) + 1L]] <-
        data.frame(
          read_id = read_id,
          type = event_type,
          start = event_start,
          end = event_end,
          position_after = NA_integer_,
          ref = ref_allele,
          alt = alt_allele,
          label = paste0(
            event_type, "_", event_start, "_", event_end,
            ":", ref_allele, ">", alt_allele
          ),

          trim_region_start = region_start,
          trim_region_end = region_end,
          reference_aligned_trim = ref_read$reference_aligned,
          read_aligned_trim = ref_read$read_aligned,

          stringsAsFactors = FALSE
        )
    }
  }

  all_events <- c(substitution_events, indel_events)

  # 没有变异时输出 WT
  if (length(all_events) == 0L) {
    event_table <- data.frame(
      read_id = read_id,
      type = "WT",
      start = NA_integer_,
      end = NA_integer_,
      position_after = NA_integer_,
      ref = NA_character_,
      alt = NA_character_,
      label = "WT",

      trim_region_start = region_start,
      trim_region_end = region_end,
      reference_aligned_trim = ref_read$reference_aligned,
      read_aligned_trim = ref_read$read_aligned,

      stringsAsFactors = FALSE
    )

    return(list(
      status = "PASS",
      read_id = read_id,
      combination = "WT",
      events = event_table
    ))
  }

  event_table <- do.call(rbind, all_events)

  # 如果同一位置同时出现插入和缺失，标记为 OTHER。
  # 原始 IN/DEL 记录仍保留在事件表中，便于检查。
  insertion_positions <- event_table$position_after[event_table$type == "IN"]
  deletion_starts <- event_table$start[event_table$type == "DEL"]

  if (length(intersect(insertion_positions, deletion_starts)) > 0L) {
    event_table$type <- "OTHER"
    event_table$label <- paste0("OTHER{", event_table$label, "}")
  }

  event_table$read_id <- read_id
  event_table <- event_table[
    , c("read_id", "type", "start", "end",
        "position_after", "ref", "alt", "label",
        "trim_region_start", "trim_region_end",
        "reference_aligned_trim", "read_aligned_trim")
  ]

  # 按参考位置排序；插入事件按其所在边界排序
  sort_pos <- ifelse(
    is.na(event_table$position_after),
    event_table$start,
    event_table$position_after + 0.5
  )
  event_table <- event_table[order(sort_pos), , drop = FALSE]

  list(
    status = "PASS",
    read_id = read_id,
    combination = paste(event_table$label, collapse = ";"),
    events = event_table
  )
}

filter_reads_by_region <- function(
    aln,
    region_start,
    region_end,
    min_anchor = 2
) {
  # 解析每条 CIGAR 中消耗 reference 坐标的部分
  ref_length <- vapply(as.character(aln$cigar), function(cigar) {
    ops <- regmatches(
      cigar,
      gregexpr("[0-9]+[MIDNSHP=X]", cigar)
    )[[1]]

    lens <- as.integer(sub("[A-Z=]$", "", ops))
    codes <- sub("^[0-9]+", "", ops)

    sum(lens[codes %in% c("M", "D", "N", "=", "X")])
  }, integer(1))

  read_start <- as.integer(aln$pos)
  read_end <- read_start + ref_length - 1L

  required_start <- max(1L, region_start - min_anchor)
  required_end <- region_end + min_anchor

  keep <- read_start <= required_start &
    read_end >= required_end

  aln[keep,,drop = FALSE]
}

# 第一步：逐 read 解析一次，生成事件表
make_variant_events <- function(
    aln_valid,
    reference = ref,
    region_start = 50,
    region_end = 220,
    min_anchor = 2,
    sample_id = "sample",
    mc.cores = NULL
) {
  reference_name <- names(reference)

  if (length(reference_name) != 1L) {
    stop("reference 必须是一个带名称的序列。")
  }

  if (is.null(mc.cores)) {
    ncores <- parallel::detectCores() - 1L
    if (is.na(ncores) || ncores < 1L) ncores <- 1L
  } else {
    ncores <- max(1L, as.integer(mc.cores))
  }

  aln_valid <- filter_reads_by_region(
    aln = aln_valid,
    region_start = region_start,
    region_end = region_end,
    min_anchor = min_anchor
  )

  events <- parallel::mclapply(
    seq_along(aln_valid$qname),
    function(i) {
      one_read <- process_one_read2(
        read_id = aln_valid$qname[i],
        read_seq = as.character(aln_valid$seq[i]),
        cigar = as.character(aln_valid$cigar[i]),
        ref_start = as.integer(aln_valid$pos[i]),
        reference = reference,
        region_start = region_start,
        region_end = region_end,
        min_anchor = min_anchor
      )

      one_read$events
    },
    mc.cores = ncores
  )

  events <- Filter(Negate(is.null), events)

  if (length(events) == 0L) {
    return(data.frame())
  }

  events_all <- do.call(rbind, events)

  # 给结果附上样本和分析区域信息，方便合并不同 BAM
  events_all$sample_id <- sample_id
  events_all$Ref_name <- reference_name
  events_all$region_start <- region_start
  events_all$region_end <- region_end

  rownames(events_all) <- NULL
  events_all
}


# 第二步：对已生成的事件表汇总，并按 min_ratio 过滤
summarise_variant_events <- function(events_all, min_ratio = 1) {
  if (min_ratio > 20) {
    warning("min_ratio = ", min_ratio, "%，请确认这个阈值是否符合预期。")
  }

  empty_summary <- data.frame(
    sample_id = character(),
    Ref_name = character(),
    Index = integer(),
    Read_count = integer(),
    Ratio = numeric(),
    Variation_count = integer(),
    Variation_combination = character(),
    region_start = integer(),
    region_end = integer(),
    reference_aligned_trim = character(),
    read_aligned_trim = character(),
    stringsAsFactors = FALSE
  )

  if (nrow(events_all) == 0L) {
    return(list(
      events_all = events_all,
      event_frequency = data.frame(),
      var_info = empty_summary,
      var_info_filtered = empty_summary
    ))
  }

  # 兼容旧事件表：没有 sample_id 时，全部视为一个样本
  if (!"sample_id" %in% names(events_all)) {
    events_all$sample_id <- "sample"
  }

  events_exclude_no_call <- events_all[
    events_all$label != "NO_CALL",
    ,
    drop = FALSE
  ]

  if (nrow(events_exclude_no_call) == 0L) {
    return(list(
      events_all = events_all,
      event_frequency = data.frame(),
      var_info = empty_summary,
      var_info_filtered = empty_summary
    ))
  }

  # 每个 label 在每个样本、参考序列中支持的 read 数
  support_data <- unique(events_exclude_no_call[
    , c("sample_id", "Ref_name", "read_id", "label"),
    drop = FALSE
  ])

  event_frequency <- stats::aggregate(
    read_id ~ sample_id + Ref_name + label,
    data = support_data,
    FUN = length
  )
  names(event_frequency)[names(event_frequency) == "read_id"] <-
    "Supporting_reads"

  callable_reads <- unique(events_exclude_no_call[
    , c("sample_id", "Ref_name", "read_id"),
    drop = FALSE
  ])

  denominators <- stats::aggregate(
    read_id ~ sample_id + Ref_name,
    data = callable_reads,
    FUN = length
  )
  names(denominators)[names(denominators) == "read_id"] <-
    "Callable_reads"

  event_frequency <- merge(
    event_frequency,
    denominators,
    by = c("sample_id", "Ref_name"),
    all.x = TRUE
  )
  event_frequency$Ratio <- 100 *
    event_frequency$Supporting_reads / event_frequency$Callable_reads

  # 低频 label 不纳入组合
  frequent_labels <- event_frequency[
    event_frequency$Ratio >= min_ratio,
    c("sample_id", "Ref_name", "label"),
    drop = FALSE
  ]

  frequent_key <- paste(
    frequent_labels$sample_id,
    frequent_labels$Ref_name,
    frequent_labels$label,
    sep = "\034"
  )

  event_key <- paste(
    events_exclude_no_call$sample_id,
    events_exclude_no_call$Ref_name,
    events_exclude_no_call$label,
    sep = "\034"
  )

  keep_event <- event_key %in% frequent_key

  # 按样本、参考序列和 read 汇总事件组合
  read_key <- paste(
    events_exclude_no_call$sample_id,
    events_exclude_no_call$Ref_name,
    events_exclude_no_call$read_id,
    sep = "\034"
  )

  read_groups <- split(seq_len(nrow(events_exclude_no_call)), read_key)

  get_first_column <- function(x, column) {
    if (column %in% names(x)) as.character(x[[column]][1]) else NA_character_
  }

  comb <- lapply(read_groups, function(idx) {
    ev <- events_exclude_no_call[idx, , drop = FALSE]
    retained_labels <- unique(ev$label[keep_event[idx]])

    if (length(retained_labels) == 0L) {
      retained_labels <- "WT"
    }

    # 按事件在表中的顺序保留组合标签
    combination <- paste(retained_labels, collapse = ";")
    variation_count <- sum(retained_labels != "WT")

    data.frame(
      sample_id = as.character(ev$sample_id[1]),
      Ref_name = as.character(ev$Ref_name[1]),
      read_id = as.character(ev$read_id[1]),
      Variation_combination = combination,
      Variation_count = variation_count,
      region_start = as.integer(ev$region_start[1]),
      region_end = as.integer(ev$region_end[1]),
      reference_aligned_trim =
        get_first_column(ev, "reference_aligned_trim"),
      read_aligned_trim =
        get_first_column(ev, "read_aligned_trim"),
      stringsAsFactors = FALSE
    )
  })

  comb <- do.call(rbind, comb)
  rownames(comb) <- NULL

  # 汇总每种 read 组合的数量和比例
  combo_key <- paste(
    comb$sample_id,
    comb$Ref_name,
    comb$Variation_combination,
    sep = "\034"
  )
  combo_groups <- split(seq_len(nrow(comb)), combo_key)

  mode_or_na <- function(x) {
    x <- x[!is.na(x) & nzchar(x)]
    if (length(x) == 0L) return(NA_character_)
    names(which.max(table(x)))
  }

  var_info <- lapply(combo_groups, function(idx) {
    z <- comb[idx, , drop = FALSE]

    denominator <- denominators$Callable_reads[
      denominators$sample_id == z$sample_id[1] &
        denominators$Ref_name == z$Ref_name[1]
    ][1]

    data.frame(
      sample_id = z$sample_id[1],
      Ref_name = z$Ref_name[1],
      Read_count = nrow(z),
      Ratio = round(100 * nrow(z) / denominator, 2),
      Variation_count = z$Variation_count[1],
      Variation_combination = z$Variation_combination[1],
      region_start = z$region_start[1],
      region_end = z$region_end[1],
      reference_aligned_trim = mode_or_na(z$reference_aligned_trim),
      read_aligned_trim = mode_or_na(z$read_aligned_trim),
      stringsAsFactors = FALSE
    )
  })

  var_info <- do.call(rbind, var_info)
  rownames(var_info) <- NULL

  order_idx <- order(
    var_info$sample_id,
    var_info$Ref_name,
    -var_info$Read_count
  )
  var_info <- var_info[order_idx, , drop = FALSE]

  group_key <- paste(var_info$sample_id, var_info$Ref_name, sep = "\034")
  var_info$Index <- ave(
    seq_len(nrow(var_info)),
    group_key,
    FUN = seq_along
  )

  var_info <- var_info[
    , c(
      "sample_id", "Ref_name", "Index", "Read_count", "Ratio",
      "Variation_count", "Variation_combination",
      "region_start", "region_end",
      "reference_aligned_trim", "read_aligned_trim"
    ),
    drop = FALSE
  ]

  var_info_filtered <- var_info[
    var_info$Ratio >= min_ratio,
    ,
    drop = FALSE
  ]

  list(
    events_all = events_all,
    event_frequency = event_frequency,
    var_info = var_info,
    var_info_filtered = var_info_filtered
  )
}

variant_table = function(aln_valid, reference = ref, region_start = 50,
                         region_end = 220,
                         min_anchor = 2,
                         min_ratio = 1){   # This parameter should put in a separate function to
  reference_name = names(reference)
  if(length(reference_name) != 1){
    stop("Number of reference must be 1, and it must have a name.")
  }

  if(min_ratio>20){
    warning("min_ratio is greater than ", min_ratio, ", which is very high.")
  }
  library(parallel)
  events = mclapply(seq_along(aln_valid$qname), function(mi){
    one_read <- process_one_read2(
      read_id = aln_valid$qname[mi],
      read_seq = as.character(aln_valid$seq[mi]),
      cigar = as.character(aln_valid$cigar[mi]),
      ref_start = as.integer(aln_valid$pos[mi]),
      reference = reference,
      region_start = region_start,
      region_end = region_end,
      min_anchor = min_anchor
    )
    one_read$events
  }, mc.cores = detectCores() - 1)

  events_all = do.call("rbind", events)

  events_exclude_no_call = subset(events_all, label != "NO_CALL")
  if(nrow(events_exclude_no_call) == 0){
    # tbl = data.frame(Var1 = character(), Freq = integer(), Ratio = double())
    return(list(events_all = events_all,
                tbl_comb = data.frame(Ref_name = reference_name,
                                      events_combination = "NO_CALL", Freq = nrow(events_all), Ratio = 100)
    ))
  }
  tbl = as.data.frame(table(events_exclude_no_call$label))
  tbl$Ratio = tbl$Freq/nrow(events_exclude_no_call) * 100
  tbl = subset(tbl, Ratio >= min_ratio)

  events_valid = events_exclude_no_call[events_exclude_no_call$label %in% tbl$Var1,,drop = FALSE]

  Variation_combination = tapply(events_valid$label, events_valid$read_id, function(x){paste(x, collapse = ";")})
  Variation_combination_n = tapply(events_valid$label, events_valid$read_id, length)
  comb = do.call("rbind", tapply(1:nrow(events_valid), events_valid$read_id,
                                 function(x){
                                   ev = events_valid[x, , drop = FALSE]
                                   data.frame(
                                     read_id = unique(ev$read_id),
                                     Variation_combination = paste(ev$label, collapse = ";"),
                                     Variation_count = length(x),
                                     reference_aligned_trim = unique(ev$reference_aligned_trim),
                                     read_aligned_trim = unique(ev$read_aligned_trim)
                                   )
                                 }))
  var_info = do.call("rbind", tapply(1:nrow(comb), comb$Variation_combination, function(x){
    comb_x = comb[x, , drop = FALSE]
    data.frame(Ref_name = reference_name,
               Variation_combination = unique(comb_x$Variation_combination),
               Variation_count = unique(comb_x$Variation_count),
               Read_count = length(comb_x$Variation_combination),
               region_start = region_start,
               region_end = region_end,
               reference_aligned_trim = names(tail(sort(table(comb_x$reference_aligned_trim)), 1)),
               read_aligned_trim = names(tail(sort(table(comb_x$read_aligned_trim)), 1)))
  }))
  var_info$Ratio = round(var_info$Read_count/sum(var_info$Read_count) * 100, digits = 2)
  var_info = var_info[order(var_info$Read_count, decreasing = TRUE),, drop = FALSE]
  var_info$Index = 1:nrow(var_info)
  var_info = var_info[, c("Ref_name", "Index", "Read_count", "Ratio",
                          "Variation_count", "Variation_combination",
                          "region_start", "region_end",
                          "reference_aligned_trim", "read_aligned_trim"), drop = FALSE]
  rownames(var_info) = 1:nrow(var_info)
  var_info_filtered = subset(var_info, Read_count/sum(Read_count)*100 >= min_ratio)
  var_info_filtered$Ratio = round(var_info_filtered$Read_count/sum(var_info_filtered$Read_count) * 100, digits = 2)

  return(list(events_all = events_all, var_info = var_info, var_info_filtered = var_info_filtered))
}


color_excel_by_barcode_ref <- function(
    data,
    file,
    barcode_col = "barcode",
    ref_col = "Ref_name",
    sheet = "Results",
    main_colors = c("#4E79A7", "#59A14F", "#E15759")
) {
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("请先安装 openxlsx：install.packages('openxlsx')")
  }
  if (length(main_colors) != 3L) {
    stop("main_colors 必须提供 3 个颜色。")
  }

  barcode <- as.character(data[[barcode_col]])
  ref_name <- as.character(data[[ref_col]])

  # 按首次出现顺序给 barcode 循环分配三个主色
  barcodes <- unique(barcode)
  barcode_colors <- setNames(
    rep(main_colors, length.out = length(barcodes)),
    barcodes
  )

  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheet)
  openxlsx::writeData(wb, sheet, data)

  for (bc in barcodes) {
    bc_rows <- which(barcode == bc)
    refs <- unique(ref_name[bc_rows])

    # 每种 Ref_name 使用同一主色系的不同明暗度
    white_mix <- if (length(refs) == 1L) {
      0.45
    } else {
      seq(0.78, 0.25, length.out = length(refs))
    }

    base_rgb <- grDevices::col2rgb(barcode_colors[[bc]])[, 1] / 255

    for (i in seq_along(refs)) {
      rows <- which(barcode == bc & ref_name == refs[i])

      # white_mix 越大，填充色越浅
      fill_rgb <- base_rgb * (1 - white_mix[i]) + white_mix[i]
      fill_color <- grDevices::rgb(
        fill_rgb[1], fill_rgb[2], fill_rgb[3],
        maxColorValue = 1
      )

      openxlsx::addStyle(
        wb,
        sheet = sheet,
        style = openxlsx::createStyle(fgFill = fill_color),
        rows = rows + 1L,  # Excel 第 1 行是表头
        cols = seq_len(ncol(data)),
        gridExpand = TRUE,
        stack = TRUE
      )
    }
  }

  openxlsx::freezePane(wb, sheet, firstRow = TRUE)
  openxlsx::saveWorkbook(wb, file, overwrite = TRUE)

  invisible(file)
}

run_editing_batch <- function(
    sample_info,
    ref_all,
    re_names,
    bam_root = "/data/project_delivery/20260904PCR/PM202609003NS_100bp_bkup_original/variant_results",
    cache_dir = "event_cache",
    min_anchor = 0,
    min_ratio = 1,
    mc.cores = NULL,
    force_rebuild = FALSE,
    colored_xlsx = NULL
) {
  # sample_info 可以是 Excel 路径，也可以是 data.frame
  if (is.character(sample_info) && length(sample_info) == 1L) {
    sample_info_0 <- as.data.frame(
      readxl::read_excel(sample_info)
    )
  } else {
    sample_info_0 <- as.data.frame(sample_info)
  }

  if (!"barcode" %in% names(sample_info_0)) {
    stop("sample_info 中缺少 barcode 列。")
  }

  # 兼容 Excel 原始列名和已经整理过的列名
  if (!"ref_name" %in% names(sample_info_0)) {
    sample_info_0$ref_name <- sample_info_0$Reference_name
  }
  if (!"region_start" %in% names(sample_info_0)) {
    sample_info_0$region_start <-
      sample_info_0$Reference_Region_Start
  }
  if (!"region_end" %in% names(sample_info_0)) {
    sample_info_0$region_end <-
      sample_info_0$Reference_Region_End
  }

  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  bam_param <- Rsamtools::ScanBamParam(
    flag = Rsamtools::scanBamFlag(
      isUnmappedQuery = FALSE,
      isSecondaryAlignment = FALSE,
      isSupplementaryAlignment = FALSE,
      isDuplicate = FALSE,
      isNotPassingQualityControls = FALSE
    ),
    what = c("qname", "rname", "pos", "cigar", "mapq", "seq", "flag")
  )

  output_rows <- list()
  per_sample_results <- list()
  out_i <- 1L

  safe_name <- function(x) {
    gsub("[^A-Za-z0-9_.-]", "_", x)
  }

  for (i in seq_len(nrow(sample_info_0))) {
    barcode <- as.character(sample_info_0$barcode[i])

    ref_names <- strsplit(
      as.character(sample_info_0$ref_name[i]),
      split = "\\|"
    )[[1]]
    ref_names <- trimws(ref_names)
    ref_names <- ref_names[nzchar(ref_names)]

    region_start <- as.integer(sample_info_0$region_start[i])
    region_end <- as.integer(sample_info_0$region_end[i])

    message("Processing ", barcode, ": ", paste(ref_names, collapse = ", "))

    cache_files <- setNames(
      file.path(
        cache_dir,
        paste0(
          safe_name(barcode), "__",
          safe_name(ref_names), "__",
          region_start, "_", region_end, "__anchor",
          min_anchor, ".events.rds"
        )
      ),
      ref_names
    )

    # 缓存有缺失时才读取该 barcode 的 BAM
    need_bam <- force_rebuild || any(!file.exists(cache_files))

    if (need_bam) {
      bam_file <- file.path(
        bam_root,
        barcode,
        "alignments",
        paste0(barcode, ".aligned.sorted.bam")
      )

      if (!file.exists(bam_file)) {
        warning("找不到 BAM，跳过：", bam_file)
        next
      }

      aln <- as.data.frame(
        Rsamtools::scanBam(bam_file, param = bam_param)[[1]]
      )

      original_rname <- as.character(aln$rname)
      mapped_rname <- unname(re_names[original_rname])
      mapped_rname[is.na(mapped_rname)] <-
        original_rname[is.na(mapped_rname)]
      aln$rname <- mapped_rname

      for (ref_name in ref_names) {
        cache_file <- cache_files[[ref_name]]

        if (!force_rebuild && file.exists(cache_file)) {
          next
        }

        if (!ref_name %in% names(ref_all)) {
          warning("ref_all 中找不到参考序列：", ref_name)
          next
        }

        aln_ref <- aln[
          !is.na(aln$rname) & aln$rname == ref_name,
          ,
          drop = FALSE
        ]

        if (nrow(aln_ref) > 0L) {
          aln_ref <- filter_reads_by_region(
            aln = aln_ref,
            region_start = region_start,
            region_end = region_end,
            min_anchor = min_anchor
          )
        }

        if (nrow(aln_ref) == 0L) {
          events_all <- data.frame()
        } else {
          events_all <- make_variant_events(
            aln_valid = aln_ref,
            reference = ref_all[ref_name],
            region_start = region_start,
            region_end = region_end,
            min_anchor = min_anchor,
            sample_id = barcode,
            mc.cores = mc.cores
          )
        }

        saveRDS(events_all, cache_file)
      }
    }

    barcode_result <- list()

    for (ref_name in ref_names) {
      cache_file <- cache_files[[ref_name]]

      if (!file.exists(cache_file)) {
        next
      }

      events_all <- readRDS(cache_file)

      if (nrow(events_all) == 0L) {
        next
      }

      summary_result <- summarise_variant_events(
        events_all = events_all,
        min_ratio = min_ratio
      )

      tab <- summary_result$var_info_filtered

      if (nrow(tab) == 0L) {
        next
      }

      tab$barcode <- barcode

      if ("Sample_HiTOM" %in% names(sample_info_0)) {
        tab$Sample_HiTOM <- as.character(sample_info_0$Sample_HiTOM[i])
      }

      barcode_result[[ref_name]] <- tab
      output_rows[[out_i]] <- tab
      out_i <- out_i + 1L
    }

    per_sample_results[[paste0(barcode, "_", i)]] <- barcode_result
  }

  if (length(output_rows) == 0L) {
    result_table <- data.frame()
  } else {
    result_table <- do.call(rbind, output_rows)
    rownames(result_table) <- NULL

    front_cols <- intersect(
      c("Sample_HiTOM", "barcode"),
      names(result_table)
    )

    result_table <- result_table[
      , c(front_cols, setdiff(names(result_table), front_cols)),
      drop = FALSE
    ]
  }

  if (!is.null(colored_xlsx)) {
    if (!exists("color_excel_by_barcode_ref", mode = "function")) {
      stop("请先定义 color_excel_by_barcode_ref()。")
    }

    color_excel_by_barcode_ref(
      data = result_table,
      file = colored_xlsx,
      barcode_col = "barcode",
      ref_col = "Ref_name",
      sheet = "Results"
    )
  }

  list(
    table = result_table,
    per_sample = per_sample_results
  )
}
