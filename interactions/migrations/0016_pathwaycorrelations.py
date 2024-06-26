# Generated by Django 3.2.15 on 2022-10-05 22:53

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('interactions', '0015_auto_20220913_1448'),
    ]

    operations = [
        migrations.CreateModel(
            name='pathwayCorrelations',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('p1a', models.CharField(choices=[('s', 'Sending'), ('r', 'Recieving')], default='s', max_length=1)),
                ('p2a', models.CharField(choices=[('s', 'Sending'), ('r', 'Recieving')], default='s', max_length=1)),
                ('correlation', models.FloatField(default=0)),
                ('pval', models.FloatField(default=0)),
                ('pathway1', models.ForeignKey(null=True, on_delete=django.db.models.deletion.SET_NULL, related_name='pathway1', to='interactions.pathway')),
                ('pathway2', models.ForeignKey(null=True, on_delete=django.db.models.deletion.SET_NULL, related_name='pathway2', to='interactions.pathway')),
            ],
        ),
    ]
