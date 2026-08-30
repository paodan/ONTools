

library(readxl)
split_sample_info = function(sampleInfo,
                             by = c("Project_ID", "Expected_Size_bp"),
                             keep = c("Sample_ID", "Project_ID", "Barcode_ID", "Sample_Type",
                                      "Expected_Size_bp", "Min_Read_Length", "Max_Read_Length",
                                      "Gel_lane", "Gel_Quality", "Clear_Band",
                                      "Analysis_Type", "Primer_F", "Primer_R"),
                             sheet = "2_metadata",
                             skip = NULL){
  if(is.character(sampleInfo) & file.exists(sampleInfo)){
    if(!is.null(sheet)){
      if(is.null(skip)){
        test = suppressMessages(readxl::read_excel(path = sampleInfo, sheet = sheet))
        skip = which(test[[1]] == "[data]")+1
      }
      sampleInfo_2 = as.data.frame(readxl::read_excel(path = sampleInfo, sheet = sheet, skip = skip))
    } else {
      if(is.null(skip))
        skip = 0
      sampleInfo_2 = read.csv(file = sampleInfo, header = TRUE, sep = ",", skip = skip)
    }

  } else {
    sampleInfo_2 = sampleInfo
  }

  sampleInfo_2[["Folders"]] = unname(apply(sampleInfo_2[, by], 1, paste0, collapse = "_"))
  readLen = ONTools::parse_size(sampleInfo_2[["Expected_Size_bp"]])
  sampleInfo_2[["Min_Read_Length"]] =


}

