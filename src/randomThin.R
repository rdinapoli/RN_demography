randomThin = function(bins)
{
  df = data.frame(RN=1:length(bins),bins=bins)
  allbins <- unique(as.character(bins))
  sitelist <- vector(mode="list", length=length(allbins))
  for (a in 1:length(allbins)){
    df1 <- df[df$bins==allbins[a],]
    if (nrow(df1)<=1){
      sitelist[[a]] <- df1$RN
    } else {sitelist[[a]] <- sample(df1$RN,1)}
  }
  index=sort(unlist(sitelist))
  return(index)
}