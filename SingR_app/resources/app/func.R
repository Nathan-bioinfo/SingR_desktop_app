load_packages = function()
{
  library(shiny)
  library(shinydashboard)
  library(shinyFiles)
  library(shinyjs)
  library(Seurat)
  library(ggplot2)
  library(EnhancedVolcano)
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)
  library(org.Rn.eg.db)
  library(tidyverse)
  library(MAST)
  library(fgsea)
  library(enrichplot)
  library(later)
  library(nichenetr)
  library(Hmisc)
  import::from(shinyalert,useShinyalert,shinyalert)
  import::from(shinybusy,add_busy_spinner)
  import::from(shinyalert,useShinyalert,shinyalert)
  import::from(shinydashboardPlus,updateBox,box)
  import::from(SeuratObject,AddMetaData,Idents)
  import::from(ggpubr,stat_compare_means,ggviolin)
  import::from(clusterProfiler,enrichGO,gseGO,gseaplot)
  import::from(dittoSeq,dittoPlot,dittoBarPlot,dittoDotPlot)
  import::from(cowplot,plot_grid)
  import::from(dplyr,full_join)
}

col_options = function(p_color_type,p_color_list,p_n_color)
{
  if("custom" %in% p_color_type)
  {
    return(p_color_list)
  }
  else if("white" %in% p_color_type)
  {
    white_vec = rep("white",p_n_color)
    return(white_vec)
  }
  else if("black" %in% p_color_type)
  {
    black_vec = rep("black",p_n_color)
    return(black_vec)
  }
  else if("default" %in% p_color_type)
  {
    return(NULL)
  }
}

n_color = function(p_input_group_by,p_input_idents,p_idents_activation,p_df)
{
  if(length(p_input_idents > 0) && p_idents_activation)
  {
    result = length(p_input_idents)
  }
  else
  {
    result = length(unique(eval(parse(text = paste0("p_df$",p_input_group_by)))))
  }
  return(result)
}

# check if a string exist in a list
exist_in = function(element,list)
{
  if(element %in% list){
    return(T)
  }else
  {
    return(F)
  }
}

# color list
colors = function()
{
  list = c("black","white","grey","lightgrey","dodgerblue","blue","darkblue","navy","lightblue","darkturquoise","cyan","red","darkred","firebrick","hotpink","pink","plum2","magenta","purple","violet","darkorchid","darkolivegreen1","green","lightgreen","lawngreen","burlywood","orange","darksalmon","chocolate","orange2","brown","palevioletred","yellow","lightgoldenrod")
  return(list)
}

# gs list creation and global value subset df storage
create_list = function(ref_df)
{
  my_list = c()
  for (i in 1:length(colnames(ref_df))) {
    # get input df real name
    name_df = deparse(substitute(ref_df))
    # store column of ref_df in name
    name = eval(parse(text = paste0(name_df,"$",colnames(ref_df[i]))))
    # convert human to mice genes
    name = convert_human_to_mouse_symbols(name, version = 1)
    # convert mouse to human
    # name = convert_mouse_to_human_symbols(name, version = 1)
    # remove NA
    name = name[!is.na(name)]
    # convert to dataframe
    name = as.data.frame(name)
    # put rownames (mices genes on a col)
    name[,2] = rownames(name)
    # remove rownames
    rownames(name) = NULL
    # named colnames
    colnames(name) = c("MGI.symbol","HGNC.symbol")
    # store df in global environment
    assign(colnames(ref_df[i]),name, envir = .GlobalEnv)
    # add value to a list
    my_list[[colnames(ref_df[i])]] = eval(parse(text = paste0(colnames(ref_df[i]),"$HGNC.symbol")))
    # confirmation message
    print(paste("column",colnames(ref_df[i]),"ok"))
  }
  return(my_list)
}

add_module_score = function(list,df)
{
  tryCatch(
    {
      for (i in 1:length(df)) 
      {
        GS<-df[i]
        for(j in 1:length(list))
        {
          meta = list[j]
          print(paste("meta = ",meta))
          input = AddModuleScore(
            eval(parse(text = meta)),
            GS,
            nbin = 24,
            ctrl = 500,
            k = FALSE,
            assay = NULL,
            name = names(GS),
            seed = 1,
            search = FALSE
          )
          print("test")
          assign(meta,input,envir = .GlobalEnv)
        }
      }
      print("module score ok")
      return(NULL)
    },
    error = function(e)
    {
      print("test error")
      print(paste("error = ",e))
      return(NULL)
    }
  )
}

regroup = function(df,p_list)
{
  tryCatch(
    {
      my_list = list()
      for (iter in 1:length(p_list)) {
        name = p_list[iter]
        input = tibble(names = names(eval(parse(text = name))$nFeature_RNA))
        for (j in 1:length(names(df))) {
          col = paste0(name,"$",names(df[j]),"1")
          input[,j+1] = eval(parse(text = col))
        }
        colnames(input) = c("cells",colnames(df))
        input = as.data.frame(input)
        rownames(input) = input[,1]
        input[,1] = NULL
        input[,length(input)+1] = name
        assign(paste0("list_",name),input,envir = .GlobalEnv)
        if(iter == 1)
        {
          my_list = eval(parse(text = paste0("list_",name)))
        }
        else
        {
          my_list = full_join(my_list,eval(parse(text = paste0("list_",name))))
        }
      }
      print("list created")
      return(my_list)
    },
    error = function(e)
    {
      return(message(e))
    }
  )
}

get_pairs = function(p_list)
{
  my_list = list()
  iter = 1
  for(i in 1:(length(p_list)-1))
  {
    for(j in (i+1):(length(p_list)))
    {
      my_list[[iter]] = c(p_list[i],p_list[j])
      iter = iter +1
    }
  }
  return(my_list)
}

get_pairs = function(p_list)
{
  my_list = c()
  iter = 1
  for(i in 1:(length(p_list)-1))
  {
    for(j in (i+1):(length(p_list)))
    {
      input = paste(p_list[i],"VS",p_list[j])
      my_list = append(my_list,input)
    }
  }
  return(my_list)
}

