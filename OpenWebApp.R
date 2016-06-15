library("tcltk")




library('scde')
load(file=('~/RNA_Seq/Cluster/plots/MCF10CA_3D_5d_random/MCF10CA_3D_5d_random_app.RData'))
show.app(app, 'pollen', browse = TRUE, port = 1468)
mywait <- function() {
  tt <- tktoplevel()
  tkpack( tkbutton(tt, text='Continue', command=function()tkdestroy(tt)),
          side='bottom')
  tkbind(tt,'<Key>', function()tkdestroy(tt) )

  tkwait.window(tt)
}
mywait()
