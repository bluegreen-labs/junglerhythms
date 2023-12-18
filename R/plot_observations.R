#' Plot raw annotations
#'
#' Visualizes the extracted annotations
#'
#' @param df annotations for a particular yearly section
#' @export
#'

plot_raw_annotations <- function(df){

  # set colour scheme for plotting
  colours = c('blue','green','black','purple','orange')

  # return empty plot if values are NA or null
  if (any(is.na(df)) | is.null(df)){
    plot(1,type='n',xaxt='n',yaxt='n',xlab='',ylab='',bty='n')
  }else{
    x = df
    x = cbind(rep(0,5),x,rep(0,5))
    plot(rep(1,337),type='n',xlab='DOY (48 week year)',ylab='',ylim=c(0,60),
         yaxt='n',
         main=sprintf('summary statistics (max. # class. %s out of 10)',
                      max(df, na.rm=T)))
    abline(v=seq(0,337,28),lwd=2,lty=2,col='grey')
    abline(v=168,lwd=2,lty=2,col='black')

    # add polygons to outline the shape of the raw classifications
    for (i in 1:5){
      s = x[i,]+(i*10)
      polygon(c(0:337),s,c(337:0),rep(i*10,337), col = colours[i])
    }

    # add y-axis labels
    text(y = seq(10, 50, by=10),
         par("usr")[1],
         labels = rev(c('flowering',
                        'fruiting',
                        'fruit drop',
                        'leaf dormancy',
                        'leaf turnover')),
         srt = 45, pos=2, xpd=TRUE)
  }
}
