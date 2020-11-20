
better.venn = function(a, b, c=NULL, d=NULL, names=NULL, ...) {
    require(limma)
	list.all <- union(a, b)
	list.all = union(list.all, c)
	list.all = union(list.all, d)
	
	list.mat = cbind(as.numeric(list.all %in% a), as.numeric(list.all %in% b))
	if (!is.null(c))
		list.mat = cbind(list.mat, as.numeric(list.all %in% c))
	if (!is.null(d))
		list.mat = cbind(list.mat, as.numeric(list.all %in% d))
	list.mat	
	list.venn <- vennCounts(list.mat)

	if (ncol(list.mat) == 4) 
		venn4(list.mat, names)
	else	
		drawVennDiagram(list.venn, names=names, ...)
}


my.doVennDiagram = function (a, b, c = NULL, d=NULL, names, ...) 
{
    require(limma)
    if (is.null(c)) {
        list.all <- union(a, b)
        list.mat <- matrix(0, nrow = length(list.all), ncol = 2)
        colnames(list.mat) <- c("list1", "list2")
        for (i in 1:length(list.all)) {
            list.mat[i, 1] <- list.all[i] %in% a
            list.mat[i, 2] <- list.all[i] %in% b
        }
        list.venn <- vennCounts(list.mat)
        drawVennDiagram(list.venn, names = names, mar = rep(1, 
            4), cex = 1.5, ...)
        ab <- intersect(which(list.mat[, 1] == 1), which(list.mat[, 
            2] == 1))
        list.ab <- vector(length = length(list.all))
        list.ab[ab] <- TRUE
        fileName <- "Venn_list_1+2.csv"
        if (!missing(names)) {
            fileName <- paste("Venn", names[1], names[2], ".csv", 
                sep = "_")
        }
        write.table(list.all[list.ab], file = fileName, sep = ",")
        print(paste("List information is written to file", fileName))
        invisible(list.all[list.ab])
    }
    else {
        list.all <- union(a, union(b, c))
        list.mat <- matrix(0, nrow = length(list.all), ncol = 3)
        colnames(list.mat) <- c("list1", "list2", "list3")
        for (i in 1:length(list.all)) {
            list.mat[i, 1] <- list.all[i] %in% a
            list.mat[i, 2] <- list.all[i] %in% b
            list.mat[i, 3] <- list.all[i] %in% c
        }
        list.venn <- vennCounts(list.mat)
        drawVennDiagram(list.venn, names = names, mar = rep(1, 
            4), cex = 1.5, ...)
        ab <- intersect(which(list.mat[, 1] == 1), which(list.mat[, 
            2] == 1))
        ac <- intersect(which(list.mat[, 1] == 1), which(list.mat[, 
            3] == 1))
        bc <- intersect(which(list.mat[, 2] == 1), which(list.mat[, 
            3] == 1))
        abc <- intersect(which(list.mat[, 1] == 1), bc)
        vd <- matrix(0, nrow = length(list.all), ncol = 4)
        rownames(vd) <- list.all
        vd[, 1][ab] <- 1
        vd[, 2][ac] <- 1
        vd[, 3][bc] <- 1
        vd[, 4][abc] <- 1
        colnames(vd) <- c("list 1_2", "list 1_3", "list 2_3", 
            "list 1_2_3")
        if (!missing(names)) {
            colnames(vd) <- c(paste(names[1], names[2], sep = "+"), 
                paste(names[1], names[3], sep = "+"), paste(names[2], 
                  names[3], sep = "+"), paste(names[1], names[2], 
                  names[3], sep = "+"))
        }
        fileName <- "Venn_list_1+2+3.csv"
        if (!missing(names)) {
            fileName <- paste("Venn", names[1], names[2], names[3], 
                ".csv", sep = "_")
        }
        write.table(vd, file = fileName, sep = ",", col.names = NA)
        print(paste("List information is written to file", fileName))
        invisible(vd)
    }
}



venn4<-function(vo,names,mar=rep(1,4)) {
	dimvo<-dim(vo)
	print(dimvo)
	if(dimvo[2] != 4)
		stop("Usage: venn4(vo,mar=rep(1,4))/n/twhere vo has 4 columns")
		plot(0,0,type="n",xlim=c(-2,2),ylim=c(-2,2),xlab="",ylab="",axes=FALSE)
		par(mar=mar)
		r = 1
		show.circle(-0.5,0.5,r)
		show.circle(0.5,0.5,r)
		show.circle(-0.5,-0.5,r)
		show.circle(0.5,-0.5,r)
		outside<-sum(apply(vo,1,sum) == 0)
		text(-1*r/1.5,1.7*r,names[1])
		text(1*r/1.5,1.7*r,names[2])
		text(-1*r/1.5,-1.7*r,names[3])
		text(1*r/1.5,-1.7*r,names[4])
#		text(0,1.4,"AB")
#		text(-1.3,0.15,"AC")
#		text(1.2,0.15,"BD")
#		text(0,-1.1,"CD")
#		text(-0.6,0.9,"ABC")
#		text(0.6,0.9,"ABD")
#		text(-0.83,-0.5,"ACD")
#		text(0.83,-0.5,"BCD")
#		text(0,0.1,"ABCD")
		par(xpd=TRUE,cex=1.5)
		text(0,-2.4,outside)
		onlyone<-apply(vo,1,sum) == 1
		A<-sum(onlyone & vo[,1])
		text(-1.2*r/1.5,1.1*r/1.5,A)
		B<-sum(onlyone & vo[,2])
		text(1.2*r/1.5,1.1*r/1.5,B)
		C<-sum(onlyone & vo[,3])
		text(-1.2*r/1.5,-1.4*r/1.5,C)
		D<-sum(onlyone & vo[,4])
		text(1.2*r/1.5,-1.4*r/1.5,D)
		onlytwo<-apply(vo,1,sum) == 2
		AB<-sum(onlytwo & vo[,1] & vo[,2])
		text(0,1.1*r/1.5,AB)
		AC<-sum(onlytwo & vo[,1] & vo[,3])
		text(-1.3*r/1.5,0*r/1.5,AC)
		BD<-sum(onlytwo & vo[,2] & vo[,4])
		text(1.2*r/1.5,-0*r/2,BD)
		CD<-sum(onlytwo & vo[,3] & vo[,4])
		text(0,-1.4*r/1.5,CD)
		onlythree<-apply(vo,1,sum) == 3
		ABC<-sum(onlythree & vo[,1] & vo[,2] & vo[,3])
		text(-0.5*r/1.5,0.5*r/1.5,ABC)
		ABD<-sum(onlythree & vo[,1] & vo[,2] & vo[,4])
		text(0.5*r/1.5,0.5*r/1.5,ABD)
		ACD<-sum(onlythree & vo[,1] & vo[,3] & vo[,4])
		text(-0.5*r/1.5,-0.5*r/1.5,ACD)
		BCD<-sum(onlythree & vo[,2] & vo[,3] & vo[,4])
		text(0.5*r/1.5,-0.5*r/1.5,BCD)
		ABCD<-sum(apply(vo,1,sum) == 4)
		text(0,0*r/1.5,ABCD)
		par(mar=c(5,4,4,2)+0.1,xpd=FALSE,cex=1)
}
drawVennDiagram = function (object, names, mar = rep(0.5, 4), cex = 1.2, cols=col.string, ...) 
{
    colors <- cols
    if (!is(object, "VennCounts")) 
        object <- vennCounts(object)
    nsets <- ncol(object) - 1
    if (nsets > 3) 
        stop("Can't plot Venn diagram for more than 3 sets")
    if (missing(names)) 
        names <- colnames(object)[1:nsets]
    counts <- object[, "Counts"]
    totalCounts = sum(counts)
    theta <- 2 * pi * (1:360)/360
    xcentres <- list(0, c(-1, 1)*1.5, c(-1, 1, 0))[[nsets]]
    ycentres <- list(0, c(0, 0), c(1/sqrt(5), 1/sqrt(5), -2/sqrt(5)))[[nsets]]
    r <- c(1.6, 1.6, 1.6)[nsets] #* 1.5 

    xtext <- list(-1.2, c(-1.2, 1.2)*1.5, c(-1.2, 1.2, 0)*1.5)[[nsets]]
    ytext <- list(1.8, c(1.8, 1.8)*1.6, c(2.4, 2.4, -3))[[nsets]]
    #old.par <- par(mar = mar)
    #on.exit(par(old.par))
	par(mar=mar)
    plot(x = 0, y = 0, type = "n", xlim = c(-4, 4), ylim = c(-4, 
        4), xlab = "", ylab = "", axes = FALSE, ...)
    for (circle in 1:nsets) {
        lines(xcentres[circle] + r * cos(theta), ycentres[circle] + 
            r * sin(theta), col = colors[circle], lwd = 3)
        text(xtext[circle], ytext[circle], names[circle], cex = cex)
    }
    switch(nsets, {
        #rect(-3, -2.5, 3, 2.5)
        text(2.3, -2.1, totalCounts, cex = cex)
        text(0, 0, counts[2], cex = cex)
    }, {
        #rect(-3, -2.5, 3, 2.5)
        text(2.3*1.6, -2.1, totalCounts, cex = cex)
        text(1.5*1.5, 0.1, counts[2], cex = cex)
        text(-1.5*1.5, 0.1, counts[3], cex = cex)
        text(0, 0.1, counts[4], cex = cex)
    }, {
        #rect(-3, -3.5, 3, 3.3)
        text(2.5, -3, totalCounts, cex = cex)
        text(0, -1.7, counts[2], cex = cex)
        text(1.5, 1, counts[3], cex = cex)
        text(0.75, -0.35, counts[4], cex = cex)
        text(-1.5, 1, counts[5], cex = cex)
        text(-0.75, -0.35, counts[6], cex = cex)
        text(0, 0.9, counts[7], cex = cex)
        text(0, 0, counts[8], cex = cex)
    })
    invisible()
}
show.circle<-function(x,y,radius,border=NULL,col=NA) {
angles<-seq(0,2*pi,by=0.04*pi)
xpos<-cos(angles)*radius+x
ypos<-sin(angles)*radius+y
polygon(xpos,ypos)
}
