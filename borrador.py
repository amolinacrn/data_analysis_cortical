            #self.CV.append(np.std(count)/np.mean(count))


            '''
            x=edges[:-1]; y=count
            x=x[y!=0]; y=y[y!=0] 
            
            if(len(y)>nN-1 and any(y)==True):
                distrib.append(y[:-1].tolist())  
                yy=y[:-1]
            elif(any(y)==True):
                distrib.append(y.tolist())
                yy=y
            else: 
                print("vector vacio")
            
            xx=np.linspace(0.05,9.95,nN-1)
            xcosim=np.zeros((len(x),len(yy)))
            cos_sim=np.zeros((len(x),len(yy)))
            for i in range(len(yy)-1):
                va=[xx[i], yy[i]]
                norm_a=np.sqrt(np.dot(va,va))
                for j in range(i+1,len(yy)):
                    vb=[xx[j], yy[j]]
                    norm_b=np.sqrt(np.dot(vb,vb))
                    pp=norm_a*norm_b
                    rest=(np.dot(va,vb)/pp)
                    cos_sim[i,j]=1-(rest) 
                    cos_sim[j,i]=1-(rest)
            '''
            # plt.plot((xx),(zx),"o")
            # xc=np.mean(cos_sim)
            # plt.yscale("log")
            # plt.xscale("log")
            # plt.show()
        
            #Hu=humbral(cos_sim)
            #io,jo=np.where(cos_sim<Hu)
            #xcosim[io,jo]=cos_sim[io,jo]
            #yi,xi=np.histogram(xcosim,bins=100)
            #plt.imshow(xcosim)
            #plt.plot(xi[:-1],yi,"o")
            #plt.yscale("log")
            #plt.xscale("log")
            
            #G_fb = nx.from_numpy_array(xcosim)
            #nx.draw_networkx(G_fb, node_size=10, with_labels = False,width=0.05)
            #plt.show()
        
            '''
            xG = nx.from_numpy_array(xcosim) #ma

            pos_p = nx.spring_layout(xG)

            part = community_louvain.best_partition(xG, weight="weight")
            resp=len(set(part.values()).union())
            nodos_por_comundad=[]
            values = list(part.values())
            nNCom=list(set(values))
                
            for cluster_id in nNCom:
                cluster = [node for node in xG.nodes() if part[node] == cluster_id]
                nodos_por_comundad.append(cluster)

            nx.draw_networkx(xG, pos=pos_p, cmap=plt.get_cmap("jet"), 
                            node_color=values, node_size=10, with_labels = False,width=0.05)
            plt.show()
            qQ=nx.community.modularity(xG,nodos_por_comundad)
            print("Modularidad = ",qQ)
            xgrade=mean_degree_network(xcosim)
            fr,x=np.histogram(xgrade.flatten(),bins=30)
            '''
            #xgrade=mean_degree_network(xcosim)
            #yi,xi=np.histogram(xgrade,bins=5)
            
            #plt.plot(xi[:-1],yi,"o");plt.show()
            #mdeg,Qq=modularidad(xcosim)
            #Mean_Degree.append(xgrade)
            #Modul.append(Qq)
            #comunid.append(mdeg)
            #meansimil.append(np.mean(xcosim))
        
        #self.CV = np.array(self.CV)
        #plt.plot(self.CV,meansimil,"o")
        #plt.yscale("log")
        #plt.xscale("log")
        #plt.show()
        '''
        fxile_txt="datos_spikes/"+direxp+"D.txt"
        fxdata = open(fxile_txt, 'w')

        try:
            for r in range(len(self.CV)):
                fxdata.write(str(self.CV[r])+"\t"+str(Mean_Degree[r])+"\t"+str(Modul[r])+"\t"+str(meansimil[r])+"\t"+str(comunid[r])+"\n")
        finally:
            fxdata.close() 
        '''
    