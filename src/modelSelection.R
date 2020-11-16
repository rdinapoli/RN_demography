modelSelection = function(x,colname='model')
{
  prop = table(x[colname])
  labels = names(prop)
  model.proportions = prop/sum(prop)
  bf.matrix=t(sapply(1:length(labels),function(x,y){y[x]/y},y=model.proportions))
  colnames(bf.matrix)=labels
  row.names(bf.matrix)=labels
  return(list(model.proportions=model.proportions,bf=bf.matrix))
}