from django.urls import path
from . import views

urlpatterns = [
    path('', views.index, name = 'index'),
    path('about/', views.help, name = 'about'),
    path('pathways/', views.PathwayListView.as_view(), name = 'pathways'),
    #path('pathways/', views.pathway_list_view, name = 'pathways'),
    path('pathways/<int:pk>', views.PathwayDetailView.as_view(), name = 'pathway-detail'),
    path('cellclasses/', views.CellClassListView.as_view(), name = 'cellclasses'),
    path('cellclasses/<int:pk>', views.CellClassDetailView.as_view(), name = 'cellclass-detail'),
    path('genes/', views.genes, name = 'genes'),
    path('ligand/<int:pk>', views.LigandDetailView.as_view(), name = 'ligand-detail'),
    path('receptor/<int:pk>', views.ReceptorDetailView.as_view(), name = 'receptor-detail'),
    path('correlations/', views.correlations, name = 'correlations'),
    path('downloadLinks', views.downloadPage, name = 'downloads'),
    path('signalingNetworks', views.signalingNetworkPage, name = 'signalingNetworklink'),
    path('correlations_heatmap/', views.correlations_heatmap, name = 'correlation_heatmap'),
    path('correlation_scatter/<pos_neg>', views.plot_correlation_scatter, name = 'correlation_scatter'),
]